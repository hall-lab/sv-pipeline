#!/bin/bash

set -eo pipefail

export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH"
RESOURCE_DIR=/opt/ccdg/cromwell/resources
#RESOURCE_DIR='.'
CROMWELL_CONF_TEMPLATE=$RESOURCE_DIR/application.conf.template
MYSQL_CONF_TEMPLATE=$RESOURCE_DIR/mysql.cnf.template

MYSQLD_PID=''
LOCK_ACQUIRED=false

MAIN_DIR="$1"

function clean_directory {
    local cleaned=$(echo "$1" | sed 's|/*$||')
    local abspath=$(cd "$cleaned" && pwd -P)
    echo "$abspath"
}

function wait_for_file {
    local file="$1"
    for i in `seq 1 28`; do
        if [[ ! -e "$file" ]]; then
            sleep $i
        fi
    done
}

function cromwell_conf {
    local dir="$1"
    echo "$dir/application.conf"
}


function set_up_conf {
    local dir="$1"
    local conf_file=$(cromwell_conf "$dir")
    if [[ -s "$conf_file" ]]; then
        echo "Using existing cromwell config at $conf_file" >&2
    else
        cat $CROMWELL_CONF_TEMPLATE | sed "s|%%SHARED_FS_DIRECTORY%%|$dir|" > "$conf_file"
        echo "Created cromwell config $conf_file" >&2
    fi
}

function has_db {
    local dir="$1"
    if [[ -d "$dir/db/run/mysqld" && -d "$dir/db/lib/mysql" && -d "$dir/db/log/mysql" && -e "$dir/db/lib/mysql/cromwell" ]]; then
        true
    else
        false
    fi
}

function create_mysql_directories {
    local dir="$1"
    for new_dir in "$dir/db/run/mysqld"  "$dir/db/lib/mysql" "$dir/db/log/mysql"
    do
        mkdir -p "$new_dir"
    done
    # TODO Are we really sure we need/want to do this?
    touch "$dir/db/log/mysql/error.log"
    chmod -R 777 "$dir/db"
}

function setup_new_database {
    echo "create database cromwell; create user 'cromwell'@'localhost' identified by 'test4cromwell'; grant all privileges on *.* to 'cromwell'@localhost;" | mysql -u root --socket=/tmp/mysqld.sock
}

function start_mysql {
    local dir="$1"
    local mysql_cnf_file=$(mysql_conf "$dir")
    mysqld_safe --defaults-file="$mysql_cnf_file" &
    MYSQLD_PID="$!"
    wait_for_file "/tmp/mysqld.sock"
}

function shutdown_mysql {
    echo "Shutting down mysql" >&2
    /usr/bin/mysqladmin -u root --socket /tmp/mysqld.sock shutdown
    MYSQLD_PID=''
}

function install_db {
    local dir="$1"
    local ldata="$dir/db/lib/mysql"
    create_mysql_directories "$dir"
    mysql_install_db --user=$USER --basedir=/usr/ --ldata=$ldata
    start_mysql "$dir"
    setup_new_database "$dir"
}

function mysql_conf {
    local dir="$1"
    echo "$dir/mysql.cnf"
}

function set_up_mysql_cnf {
    local dir="$1"
    local conf=$(mysql_conf "$dir")
    if [[ -s "$conf" ]]; then
        echo "Using existing mysql config at $conf" >&2
    else
        cat $MYSQL_CONF_TEMPLATE | sed "s|%%SHARED_FS_DIRECTORY%%|$dir|" > "$conf"
        echo "Created mysql config $conf" >&2
    fi
}

function is_locked {
    local dir="$1"
    if [[ -d "$dir/.lock" ]]; then
        true
    else
        false
    fi
}

function lock {
    local dir="$1"
    if mkdir "$dir/.lock"; then
        LOCK_ACQUIRED=true
        echo "Locked $dir" >&2
    else
        echo "Unable to lock $dir" >&2
        exit 1
    fi
}

function unlock {
    local dir="$1"
    if rmdir "$dir/.lock"; then
        LOCK_ACQUIRED=false
        echo "Unlocked $dir" >&2
    else
        echo "Unable to unlock $dir" >&2
        exit 1
    fi
}

function run_cromwell {
    local dir="$1"
    local cromwell_conf=$(cromwell_conf "$dir")
    /usr/bin/java -Dconfig.file="$cromwell_conf" -jar /app/cromwell.jar run "${@:2}"
}

function cleanup {
    local dir="$1"
    if [[ $MYSQLD_PID ]]; then
        shutdown_mysql
    fi
    if $LOCK_ACQUIRED; then
        unlock "$dir"
    fi
}

function main {
   local dir="$1"
   if [[ -d "$dir" ]]
    then
        local clean_dir=$(clean_directory "$dir")
        trap 'cleanup $(clean_directory "$MAIN_DIR")' EXIT SIGTERM SIGINT
        lock "$clean_dir"
        set_up_conf "$clean_dir"
        set_up_mysql_cnf "$clean_dir"
        if ! $(has_db "$clean_dir"); then
            # note that this also starts mysqld
            install_db "$clean_dir"
        else
            start_mysql "$clean_dir"
        fi
        echo "${@:2}" >&2
        run_cromwell "$clean_dir" "${@:2}"
    else
        echo "$dir is not a directory" >&2
        exit 1
    fi
}

main "${@}";

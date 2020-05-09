# Encrypting / Decrypting `cromwell.sv.yaml`

The `cromwell.sv.yaml` is a [YAML][1] file is used by the [cromwell-deployment][0] system.  The YAML file contains configuration information for setting up a google compute instance, MySQL database, and the software version of cromwell to install.  The file generally looks like so:

    imports:
      - path: cromwell.jinja
    
    resources:
      - name: cromwell
        type: cromwell.jinja
        properties:
          # Commented out properties have defaults set in the schema
          region:                          us-central1
          zone:                            us-central1-c
          service_account_email:           <your-service-account>
          cromwell_version:                50
          cromwell_server_machine_type:    n1-highmem-16 # 104G
          #cromwell_server_boot_disk_size:  10
          ssh_source_ranges:               [ <your-ip-range-in-cidr-format> ]
    
          cromwell_cloudsql_instance_type: db-n1-standard-16 # 60G
          cromwell_cloudsql_initial_size:  100
          cromwell_database_name:          <your-mysql-database-name>
          cromwell_cloudsql_password:      <your-mysql-database-password>
    
          labels:
            user:     <your-user-name>
            project:  <your-project-name>
            pipeline: <your-pipeline-name>

Unfortunately, the `cromwell.sv.yaml` file contains some sensitive information like the database name/password, source IP ranges, and service account information.  Because this configuration file may be placed in a public repository, we've encrypted this file with symmetric encryption via [gpg][2] as `cromwell.sv.yaml.gpg`.

To decrypt this file on your local machine, please do the following:

    cd config/cromwell-server
    gpg --output cromwell.sv.yaml --no-symkey-cache --decrypt cromwell.sv.yaml.gpg

You will be prompted for a passphrase.  Please ask prior contributors of this repository for that passphrase.  Afterwards, you should see the unencrypted `cromwell.sv.yaml` file in the directory. Now you can go about using the [cromwell-deployment][0] system to start up the server.

**NOTE:  Please do not commit and push the unencrypted `cromwell.sv.yaml` file to the public git repository!!!**

Should you need to alter the YAML file and re-encrypt it.  Please do the following:

    gpg --symmetric --cipher-algo AES256 --no-symkey-cache cromwell.sv.yaml

You will be prompted for a passphrase.  Please use the same passphrase as the one you used to originally decrypt the file.

[0]: https://github.com/hall-lab/cromwell-deployment
[1]: https://yaml.org/
[2]: https://gnupg.org/

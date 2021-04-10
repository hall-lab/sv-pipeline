#!/bin/bash

ID=$1

curl --silent -X GET --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/${ID}/status" | python -m json.tool

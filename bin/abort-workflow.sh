#!/bin/bash

ID=$1
curl -X POST --header "Content-Type: application/json" --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/${ID}/abort"

#!/bin/bash
echo "Starting..."
sed -i -e 's/github.com\/ThacherObservatory/gitlab.com\/thacher-astronomy/g' `find . -name "config"`
echo "Complete!"
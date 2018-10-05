#!/bin/bash

git config --global user.email "ypar@pennmedicine.upenn.edu"
git config --global user.name "ypar"

git checkout master
git add -u
git commit -m "update data tables: Travis Build $TRAVIS_BUILD_NUMBER [skip ci]"

git remote add deploy https://${GITHUB_TOKEN}@github.com/ypar/cimr.git > /dev/null 2>&1
git push --quiet deploy master > /dev/null 2>&1

# Create a new release to trigger Zenodo archiving
version=$(cat version.txt)
git tag $version
git push --quiet deploy --tags > /dev/null 2>&1
curl -v -i -X POST -H "Content-Type:application/json" -H "Authorization: token $GITHUB_TOKEN" https://api.github.com/repos/ypar/cimr/releases -d "{\"tag_name\":\"$version\"}"


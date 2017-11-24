#! /bin/sh

git pull github master
git add .
git commit -m "$(date)"
git push github master

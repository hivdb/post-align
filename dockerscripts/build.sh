#! /bin/bash

set -e

cp -r /app/post-align /app/local-post-align
pushd /app/local-post-align
pushd postalign
find . -name "*.c" -o -name "*.so" | xargs rm
popd
rm -rf build dist
python3.9 setup.py build_ext --inplace
find postalign -name "*.so" | while read file; do
  rm ${file%.cpython*.so}.py
done
pyinstaller postalign/entry.py -n postalign
rm -rf /app/post-align/dist/linux-amd64
mkdir -p /app/post-align/dist/
cp -r /app/local-post-align/dist /app/post-align/dist/linux-amd64

#!/bin/bash
cd ..
make wasm
cd web
cp ../nmie.js ./
cp ../nmie.wasm ./
emrun index.html

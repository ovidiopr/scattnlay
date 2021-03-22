#!/bin/bash
#cp ../nmiejs.js src
#cp ../nmiejs.wasm public
npm install
npm audit fix
npm run serve

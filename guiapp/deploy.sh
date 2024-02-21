#!/bin/bash
# build, copy build files, and copy (manually) content of index.html to https://physics.itmo.ru/ru/node/2611/edit

quasar build
# cp dist/spa/nmiejs.wasm dist/spa/js/nmiejs.wasm

cd ~/coding/kostyfisik.github.io
git reset --hard c66101abdd5db0a4ff5ef69c3328cdfdbadaa398
cp -r ../scattnlay/guiapp/dist/spa/* ./
git add *
git add  */*
git commit -am 'update'
git push -f

# # physics.itmo.ru
# cd dist/spa/css
# for file in `ls *.css`; do sed -i 's/[.]block{display:block\!important}//g' $file; done
# cd ..
# #sed -i 's=[<]base href[=]//themes//custom//physics//mie-next// [>]==g' index.html
# cd ../..

# rsync -aue ssh --progress  dist/spa/ physics@physics.ifmo.ru:/var/www/html/physicsifmoru/web/themes/custom/physics/mie
# echo
# cat dist/spa/index.html
# echo

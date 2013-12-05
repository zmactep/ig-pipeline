rm -rf ./dist
mkdir dist

echo "Building ig-config..."
cd ./ig-config
./build_package.sh
cp ig-config-0.1-1-x86_64.pkg.tar.xz ../dist
rm -rf ig-config-0.1-1-x86_64.pkg.tar.xz ig-config.tar.gz pkg src

echo "Done. Building ig-tools..."
cd ../ig-tools
./build_package.sh
cp ig-tools-0.1-1-x86_64.pkg.tar.xz ../dist
rm -rf ./src ./pkg ig-tools-0.1-1-x86_64.pkg.tar.xz ig-tools.tar.gz

echo "Done. Building ig-backend..."
cd ../ig-backend
./build_package.sh
cp ig-backend-0.1-1-x86_64.pkg.tar.xz ../dist
rm -rf ./src ./pkg ig-backend-0.1-1-x86_64.pkg.tar.xz ig-backend.tar.gz

echo "Done. Building ig-frontend..."
cd ../ig-frontend
./build_package.sh
cp ig-frontend-0.1-1-x86_64.pkg.tar.xz ../dist
rm -rf ./src ./pkg ig-frontend-0.1-1-x86_64.pkg.tar.xz ig-frontend.tar.gz

txtred='\e[0;31m' 
txtrst='\e[0m'

echo -e "${txtred}Done. All your packages are in dist folder. Now run 'su -c ./install_all.sh' to install server${txtrst}"

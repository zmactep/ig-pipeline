echo "Installing packages..."

pacman -U dist/ig-config-0.1-1-x86_64.pkg.tar.xz dist/ig-tools-0.1-1-x86_64.pkg.tar.xz dist/ig-backend-0.1-1-x86_64.pkg.tar.xz dist/ig-frontend-0.1-1-x86_64.pkg.tar.xz
echo "Done. Restarting services..."
systemctl daemon-reload
systemctl start ig-backend.service
systemctl start ig-frontend.service

echo "Done."

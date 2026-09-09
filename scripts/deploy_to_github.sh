#!/usr/bin/env bash
# Helper script to link your personal GitHub repository and push for GitHub Pages hosting

set -e

echo "=================================================================="
echo "  Deploy Navier-Stokes Laboratory to Your GitHub & GitHub Pages"
echo "=================================================================="
echo ""

# Check current git remote
CURRENT_REMOTE=$(git remote get-url origin 2>/dev/null || echo "")

echo "Current 'origin' remote is: $CURRENT_REMOTE"
echo ""
echo "To push this repository to your own GitHub account:"
echo "1. Go to https://github.com/new and create a new repository (e.g. 'navier-stokes-blowup')"
echo "2. Run the following commands:"
echo ""
echo "   # Step A: Rename upstream origin if you want to keep OpenAI's repo as upstream"
echo "   git remote rename origin upstream 2>/dev/null || true"
echo ""
echo "   # Step B: Add your GitHub repository as origin"
echo "   git remote add origin https://github.com/<YOUR_GITHUB_USERNAME>/<YOUR_REPO_NAME>.git"
echo ""
echo "   # Step C: Stage all new files and commit"
echo "   git add -A"
echo "   git commit -m 'feat: Add interactive 3D Navier-Stokes blowup simulator and docs'"
echo ""
echo "   # Step D: Push to your GitHub repository"
echo "   git push -u origin main"
echo ""
echo "=================================================================="
echo "  Enabling GitHub Pages (Free Web Hosting for LinkedIn Sharing)"
echo "=================================================================="
echo "1. Go to your repository on GitHub -> Settings -> Pages"
echo "2. Under 'Build and deployment', set Source to 'Deploy from a branch'"
echo "3. Branch: 'main', Folder: '/docs'"
echo "4. Click Save!"
echo ""
echo "Your website will be live at:"
echo "  https://<YOUR_GITHUB_USERNAME>.github.io/<YOUR_REPO_NAME>/"
echo "=================================================================="

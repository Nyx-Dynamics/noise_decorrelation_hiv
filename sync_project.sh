#!/bin/bash

# Simple sync script for noise_decorrelation_HIV project

echo "🔄 Starting Sync Process..."

# 1. Pull latest changes
echo "📥 Pulling from GitHub..."
git pull --rebase origin main

if [ $? -ne 0 ]; then
    echo "❌ Error during pull. Please resolve conflicts manually."
    exit 1
fi

# 2. Check for local changes
if [[ -n $(git status -s) ]]; then
    echo "📝 Local changes detected."
    git status -s
    
    read -p "Do you want to commit and push these changes? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        read -p "Enter commit message: " msg
        if [ -z "$msg" ]; then
            msg="Update: $(date +'%Y-%m-%d %H:%M:%S')"
        fi
        
        git add .
        git commit -m "$msg"
        
        echo "📤 Pushing to GitHub..."
        git push origin main
    fi
else
    echo "✅ No local changes to push."
fi

echo "✨ Sync Complete!"

#!/bin/bash
#
# This script creates the 'tcr_env' and then clones, patches,
# and installs tcrdist3 to avoid the "File name too long" error.
#

ENV_NAME="tcr_env"
REPO_URL="https://github.com/kmayerb/tcrdist3.git"
CLONE_DIR="./temp_tcrdist3"

echo "--- 1. Creating base environment '$ENV_NAME' from tcr_env.yml ---"
mamba env create -f tcr_env.yml

if [ $? -ne 0 ]; then
    echo "Conda environment creation failed. Aborting."
    exit 1
fi

echo "--- 2. Cloning tcrdist3 using Sparse Checkout (to skip problem file) ---"
# Clean up any previous failed attempts
rm -rf $CLONE_DIR 
mkdir $CLONE_DIR
cd $CLONE_DIR

git init -b main
git remote add origin $REPO_URL
git config core.sparseCheckout true

# Define which files/folders to check out
# This tells git:
# 1. Get all files in the root (like setup.py, MANIFEST.in, requirements.txt)
# 2. Get the main tcrdist source folder
# 3. BUT explicitly exclude the problem directory
echo "/*" > .git/info/sparse-checkout
echo "/tcrdist/*" >> .git/info/sparse-checkout
echo "!/tcrdist/data/covid19/" >> .git/info/sparse-checkout

echo "Pulling from remote... (This will skip the problem file)"
git pull origin master

# Go back to the root project directory
cd .. 

echo "--- 3. Installing tcrdist3 into '$ENV_NAME' using pip ---"
# We can now install directly. No need to patch MANIFEST.in
# because the problem folder doesn't exist on our disk.
mamba run -n $ENV_NAME pip install $CLONE_DIR

if [ $? -ne 0 ]; then
    echo "pip install failed. Aborting."
    exit 1
fi

echo "--- 4. Cleaning up temporary clone ---"
rm -rf $CLONE_DIR

echo "---"
echo "--- SUCCESS! ---"
echo "Environment '$ENV_NAME' is created and tcrdist3 is installed."
echo "To activate, run: conda activate $ENV_NAME"


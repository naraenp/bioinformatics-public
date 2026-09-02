# Instructions for environment creation

First, ensure that you have mamba and conda installed, alongside the R/Python binaries required for the subsequent scripts. Second, run all commands from the root of the repository.

```bash
# To create the relevant environments
mamba env create -f envs/preleuk_r.yml
conda activate preleuk_r

# After this check that our versions are shown properly
R --version     # should be 4.3.3
quarto check    # checks quarto to ensure it see's R/knitr correctly

# After that, please go to .vscode/settings.json.example and set the USER as per your machine and reload the window
# Ensure that VS Code root is set to the current repository, and that the R terminal itself is restarted within

# The Python environment (stages 06-07) is separate:
mamba env create -f envs/preleuk_py.yml
```

## Known issue: conda quarto package layout

The conda-forge quarto build puts its helper binaries (deno, pandoc, sass, ...) in
`<env>/bin/` while quarto expects them under `<env>/bin/tools/x86_64/`, and it
finds its resource directory one level too high. Symptoms: "No such file or
directory" for deno or pandoc, then "Module not found" for yaml-intelligence.
Fix (redo after recreating the env):

```bash
ENV=~/miniforge3/envs/preleuk_r
cd $ENV/bin/tools/x86_64
for b in deno pandoc pandoc-lua pandoc-server typst sass esbuild; do ln -sf ../../$b $b; done
mkdir -p dart-sass && ln -sf ../../../sass dart-sass/sass
mkdir -p $ENV/etc/conda/activate.d
echo "export QUARTO_SHARE_PATH=$ENV/share/quarto" > $ENV/etc/conda/activate.d/quarto-fix.sh
# re-activate the env, then: quarto check
```

For the Python stages, point quarto at the right interpreter:
`export QUARTO_PYTHON=~/miniforge3/envs/preleuk_py/bin/python`.

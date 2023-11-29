# Warp

Used to produce waved svgs of flags for Noto Emoji

Usage:

```shell
# clone noto emoji, say to ./noto-emoji
# clone warp, say to ./warp
pushd warp
# establish a python virtual environment
pip install -r requirements.txt
python wave_ninja.py --noto_dir ../noto-emoji

# your waved flags are in build/waved/
```

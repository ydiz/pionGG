#!/bin/bash
echo "Generating tags for files under this directory."
shopt -s nullglob # prevent errors caused by */*.cc not existing
ctags --c++-kinds=+p --fields=+iaS --extra=+q --language-force=C++ -f .tags -R {*,*/*}.{h,cc}

grid_path=/home/yidizhao/cooley/install/Grid/include
if [ ! -f "$HOME/cooley/.tags" ]; then
	echo "Generating tags for Grid library; output file is ~/cooley/.tags"
	ctags --c++-kinds=+p --fields=+iaS --extra=+q --language-force=C++ -f $HOME/cooley/.tags -a -R $grid_path
fi

#!/usr/bin/env bash
shopt -s globstar
clang-format -i -style=file include/**/*.h


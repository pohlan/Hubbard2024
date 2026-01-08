#!/bin/bash

ffmpeg -framerate 30 -pattern_type glob -i './frames/*.png' -c:v libx264 -r 30 output.mp4

cd "%~dp0"
ffmpeg -f image2 -r 15 -i "progress_%%5d.png" -c:v libx264 -crf 27.9 -pix_fmt yuv420p -y video.mp4

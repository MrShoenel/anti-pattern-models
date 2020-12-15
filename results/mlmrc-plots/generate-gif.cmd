cd "%~dp0"
ffmpeg -i video.mp4 -pix_fmt rgb24 -loop 0 -filter:v "scale=800:-2:flags=lanczos" -y video.gif

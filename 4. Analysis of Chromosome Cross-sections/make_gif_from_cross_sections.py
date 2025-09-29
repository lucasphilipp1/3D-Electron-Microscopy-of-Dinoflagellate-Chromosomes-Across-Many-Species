#create a gif of the XY, YZ, and XZ cross-sections of shapes along PC1 or PC2

from PIL import Image, ImageOps, ImageDraw, ImageFont
import os
import numpy as np

xy_folder = '/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC1/xy/'
yz_folder = '/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC1/yz/'
xz_folder = '/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC1/xz/'
output_gif = 'PC1_cross_section.gif'
frame_duration = 100  # ms per frame

def load_grayscale_images(folder):
    files = sorted(
        [f for f in os.listdir(folder) if f.endswith('.png')],
        key=lambda x: int(os.path.splitext(x)[0])
    )
    return [Image.open(os.path.join(folder, f)).convert("L") for f in files]

xy_imgs = load_grayscale_images(xy_folder)
yz_imgs = load_grayscale_images(yz_folder)
xz_imgs = load_grayscale_images(xz_folder)

num_frames = min(len(xy_imgs), len(yz_imgs), len(xz_imgs))
xy_imgs, yz_imgs, xz_imgs = xy_imgs[:num_frames], yz_imgs[:num_frames], xz_imgs[:num_frames]

# Pad all images to max size
def get_max_size(images):
    return max(im.width for im in images), max(im.height for im in images)

xy_max = get_max_size(xy_imgs)
xz_max = get_max_size(xz_imgs)
yz_max = get_max_size(yz_imgs)

def pad_image(im, target_size):
    return ImageOps.pad(im, target_size, method=Image.BICUBIC, color=0, centering=(0.5, 0.5))

xy_imgs = [pad_image(im, xy_max) for im in xy_imgs]
xz_imgs = [pad_image(im, xz_max) for im in xz_imgs]
yz_imgs = [pad_image(im, yz_max) for im in yz_imgs]

def get_y_centroid(im):
    arr = np.array(im)
    rows = np.where(arr > 0)
    return np.mean(rows[0]) if len(rows[0]) else im.height // 2

def align_by_centroid(im, canvas_height):
    centroid_y = get_y_centroid(im)
    shift = int(canvas_height // 2 - centroid_y)

    canvas = Image.new("L", (im.width, canvas_height), color=0)
    paste_y = max(0, shift)
    crop_y = -min(0, shift)

    region = im.crop((0, crop_y, im.width, crop_y + im.height))
    canvas.paste(region, (0, paste_y))
    return canvas

frame_height = max(xy_max[1], xz_max[1], yz_max[1]) * 2
frame_width = xy_max[0] + xz_max[0] + yz_max[0]

pc1_values = np.linspace(-20, 45, num=num_frames)

try:
    font = ImageFont.truetype("Arial.ttf", size=60)
except:
    font = ImageFont.load_default()

#generate frames of gif
frames = []
for i in range(num_frames):
    xy_aligned = align_by_centroid(xy_imgs[i], frame_height)
    xz_aligned = align_by_centroid(xz_imgs[i], frame_height)
    yz_aligned = align_by_centroid(yz_imgs[i], frame_height)

    # Create composite frame: XY | XZ | YZ
    frame = Image.new("L", (frame_width, frame_height), color=0)
    frame.paste(xy_aligned, (0, 0))
    frame.paste(xz_aligned, (xy_max[0], 0))
    frame.paste(yz_aligned, (xy_max[0] + xz_max[0], 0))

    # Add PC label
    draw = ImageDraw.Draw(frame)
    label = f"PC1: {pc1_values[i]:.1f}" #if PC1
    #label = f"PC2: {pc2_values[i]:.1f}" #if PC2
    draw.text((10, 10), label, fill=255, font=font)

    frames.append(frame)

#save gif
frames[0].save(
    output_gif,
    save_all=True,
    append_images=frames[1:],
    duration=frame_duration,
    loop=0
)

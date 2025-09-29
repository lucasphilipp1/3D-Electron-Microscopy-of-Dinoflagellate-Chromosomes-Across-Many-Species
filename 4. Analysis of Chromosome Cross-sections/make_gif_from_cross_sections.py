from PIL import Image, ImageOps, ImageDraw, ImageFont
import os
import numpy as np

xy_folder = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Spherical Harmonics Expansion/all/PC1/xy'
yz_folder = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Spherical Harmonics Expansion/all/PC1/yz/'
xz_folder = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Spherical Harmonics Expansion/all/PC1/xz/'
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

label_height = 40  # Space for XY/XZ/YZ text labels at bottom
aligned_height = max(xy_max[1], xz_max[1], yz_max[1]) * 2
frame_height = aligned_height + label_height
frame_width = xy_max[0] + xz_max[0] + yz_max[0]

try:
    font = ImageFont.truetype("Arial.ttf", size=30)
except:
    font = ImageFont.load_default()

pc1_values = np.linspace(-50, 75, num=num_frames)
pc2_values = np.linspace(-20, 45, num=num_frames)

frames = []
for i in range(num_frames):
    xy_aligned = align_by_centroid(xy_imgs[i], aligned_height)
    xz_aligned = align_by_centroid(xz_imgs[i], aligned_height)
    yz_aligned = align_by_centroid(yz_imgs[i], aligned_height)

    # Create composite frame
    frame = Image.new("L", (frame_width, frame_height), color=0)
    frame.paste(xy_aligned, (0, 0))
    frame.paste(xz_aligned, (xy_max[0], 0))
    frame.paste(yz_aligned, (xy_max[0] + xz_max[0], 0))

    draw = ImageDraw.Draw(frame)

    # Draw PC1 label at top-left
    label = f"PC1: {pc1_values[i]:.1f}"
    #label = f"PC2: {pc2_values[i]:.1f}"

    draw.text((10, 10), label, fill=255, font=font)
    
    # Draw XY, XZ, YZ labels at bottom
    y_text_pos = frame_height - label_height + 5

    # Get centers of each section
    xy_center = xy_max[0] // 2
    xz_center = xy_max[0] + xz_max[0] // 2
    yz_center = xy_max[0] + xz_max[0] + yz_max[0] // 2

    # Draw centered text (approximate with fixed offset)
    draw.text((xy_center - 15, y_text_pos), "XY", fill=255, font=font)
    draw.text((xz_center - 15, y_text_pos), "XZ", fill=255, font=font)
    draw.text((yz_center - 15, y_text_pos), "YZ", fill=255, font=font)

    frames.append(frame)

#save gif
frames[0].save(
    output_gif,
    save_all=True,
    append_images=frames[1:],
    duration=frame_duration,
    loop=0
)

file2read = ['TestSet/mask_goal_image.png'];
mask = imread(file2read);
se1 = strel('diamond', 22);
mask2 = imdilate(mask, se1);

figure, imshow(mask)
figure, imshow(mask2)

imwrite(mask2, './TestSet/mask_goal_image.bmp');
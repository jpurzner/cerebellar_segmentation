patternLength = 300; % Total length of the pattern
patternHeight = 10; % Arbitrary height to make the pattern 2D
initialPattern = zeros(patternHeight, patternLength);
initialPattern(:, 51:100) = 1;
initialPattern(:, 200:250) = 1;
rotatedPattern45 = imrotate(initialPattern, 45);
rotatedPattern90 = imrotate(initialPattern, 90);


file = '/Users/jpurzner/Dropbox/images/edu_repeat/p27/2018_05_22_s3_3_p27-0005/2018_05_22_s3_3_p27-0005_fused_crop.tif';
p27 = imread(file, 1);
p27 = mat2gray(p27);
p27 = adapthisteq(p27,'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
p27_not = not(p27 > 0.5);

% Load and normalize the dapi channel image
file = '/Users/jpurzner/Dropbox/images/edu_repeat/p27/2018_05_22_s3_3_p27-0005/2018_05_22_s3_3_p27-0005_fused_crop.tif';
dapi = imread(file, 3);
dapi = mat2gray(dapi);
dapi = adapthisteq(dapi,'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
dapi_not = not(dapi > 0.5);

neun = imread(file, 2);
neun = mat2gray(neun);
neun = adapthisteq(neun,'clipLimit',0.01,'Distribution','rayleigh', 'NumTiles', [50 50], 'Range', 'original');
neun_not = not(neun > 0.5);


filteredImage0 = conv2(double(p27), double(initialPattern), 'same');
filteredImage45 = conv2(double(p27), double(rotatedPattern45), 'same');
filteredImage90 = conv2(double(p27), double(rotatedPattern90), 'same');

threshold = 0.6; 

detectedPattern45 = filteredImage45 > threshold * max(filteredImage45(:));
detectedPattern0 = filteredImage0 > threshold * max(filteredImage45(:));
detectedPattern90 = filteredImage90 > threshold * max(filteredImage45(:));

%detectedPattern0_clean =  detectedPattern0 .* p27_not .* dapi_not .* neun_not ;
%imshow(detectedPattern0_clean)

all_p27_pattern = detectedPattern0 | detectedPattern45 | detectedPattern90;


figure, imshow(detectedPattern45 .* p27_not .* dapi_not, []), title('Detected with 45 Degrees Pattern');
figure, imshow(detectedPattern0 .* p27_not .* dapi_not, []), title('Detected with 0 Degrees Pattern');
figure, imshow(detectedPattern90 .* p27_not  .* dapi_not, []), title('Detected with 90 Degrees Pattern');

% Pattern specifics for dapi adjusted: 100 at 1, 20 at 0, 100 at 1
patternLengthDapi = 220; % Total length of the new pattern
patternHeightDapi = 10; % Height to make the pattern 2D

% Create the initial pattern for dapi
initialPatternDapi = zeros(patternHeightDapi, patternLengthDapi);
initialPatternDapi(:, 1:100) = 1; % First 100 at 1
initialPatternDapi(:, 121:220) = 1; % Second 100 at 1, after a gap of 20

% Rotate the dapi pattern by 45 and 90 degrees
rotatedPatternDapi45 = imrotate(initialPatternDapi, 45);
rotatedPatternDapi90 = imrotate(initialPatternDapi, 90);





% Apply the initial and rotated patterns as kernels to the dapi image
filteredImageDapi0 = conv2(double(dapi), double(initialPatternDapi), 'same');
filteredImageDapi45 = conv2(double(dapi), double(rotatedPatternDapi45), 'same');
filteredImageDapi90 = conv2(double(dapi), double(rotatedPatternDapi90), 'same');

% Threshold to detect the pattern
threshold = 0.5; % Adjust threshold based on the specific needs
detectedPatternDapi0 = filteredImageDapi0 > threshold * max(filteredImageDapi0(:));
detectedPatternDapi45 = filteredImageDapi45 > threshold * max(filteredImageDapi45(:));
detectedPatternDapi90 = filteredImageDapi90 > threshold * max(filteredImageDapi90(:));

% Display the results for 0, 45, and 90 degrees rotations
%figure, imshow(detectedPatternDapi0, []), title('Dapi Detected with 0 Degrees Pattern');
%figure, imshow(detectedPatternDapi45, []), title('Dapi Detected with 45 Degrees Pattern');
%figure, imshow(detectedPatternDapi90, []), title('Dapi Detected with 90 Degrees Pattern');

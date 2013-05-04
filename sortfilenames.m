function [sorted_filenames] = sortfilenames (filenames);
dates = zeros(size(filenames));
for i=1:size(filenames,1)
    f = filenames(i);
    display(f)
    t = regexpi(f,'^(?<prefix>.+)-(?<frame>\d+)\.jpg$', 'names');
    display(t{1}.frame)
    frames(i) = str2num(t{1}.frame);
end
[b,idx] = sort(frames);
sorted_filenames = filenames(idx);
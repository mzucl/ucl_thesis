function add_author_to_all_m_files()
    % Define author block
    authorBlock = sprintf([ ...
        '%% Author: Mediha Zukic\n' ...
        '%% Contact: \n' ...
        '%% Date: %s\n\n'], datestr(now, 'yyyy-mm-dd'));

    % Get all .m files in current folder
    files = dir('*.m');

    for k = 1:length(files)
        fileName = files(k).name;

        % Read current file contents
        fid = fopen(fileName, 'r');
        fileContents = fread(fid, '*char')';
        fclose(fid);

        % Check if author line already exists to avoid duplication
        if contains(fileContents, 'Author: Mediha Zukić')
            fprintf('Skipping %s (author already present)\n', fileName);
            continue;
        end

        % Prepend author block
        newContents = [authorBlock, fileContents];

        % Overwrite file with new content
        fid = fopen(fileName, 'w');
        fwrite(fid, newContents);
        fclose(fid);

        fprintf('Updated: %s\n', fileName);
    end
end
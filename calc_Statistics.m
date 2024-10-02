subdir = '';
order = 1;

threads = [1024, 512, 256, 128, 64, 32, 16, 8];
blocks=[16, 32, 64, 128, 256, 512, 1024, 2048];
		stat = [];
		for t = 1:length(blocks)
			fileName = sprintf("./%s/proc_thread_%d_block_%d_order_%d.dat", subdir, threads(t), blocks(t), order);
			fid = fopen(fileName,"r");
			if fid == -1
				fprintf('Error: Could not open file %s\n', fileName);
				continue; % Move to the next iteration of the loop
			end
			total_seconds = [];
			while(!feof(fid))
				str = fgets (fid);
				res = textscan (str, "%dm%d,%ds"); % If the 'time' command writes out the results with a decimal point instead of a decimal comma, this line should be modified to: res = textscan (str, "%dm%f");
				minutes = double(res{1});
				seconds = double(res{2});
				fraction = double(res{3}); % If the 'time' command writes out the results with a decimal point instead of a decimal comma, the decimal variable should be deleted or commented out.
				total_seconds = [total_seconds; minutes * 60 + seconds + fraction/1000 ]; % If the 'time' command writes out the results with a decimal point instead of a decimal comma, this line should be modified to: total_seconds = [total_seconds; minutes * 60 + seconds ];
			end
			fclose(fid);
			q = quantile(total_seconds, [0, 0.25, 0.5, 0.75, 1]);
			stat = [stat; blocks(t), threads(t), mean(total_seconds), q];
		end
		file_save = sprintf('grid_%d_order_%d.txt', threads(1) * blocks(1), order);
		save(file_save, 'stat', '-ascii');


threads = [1024, 512, 256, 128, 64, 32, 16, 8];
blocks=[8, 16, 32, 64, 128, 256, 512, 1024];
		stat = [];
		for t = 1:length(blocks)
			fileName = sprintf("./%s/proc_thread_%d_block_%d_order_%d.dat", subdir, threads(t), blocks(t), order);
			fid = fopen(fileName,"r");
			if fid == -1
				fprintf('Error: Could not open file %s\n', fileName);
				continue; % Move to the next iteration of the loop
			end
			total_seconds = [];
			while(!feof(fid))
				str = fgets (fid);
				res = textscan (str, "%dm%d,%ds"); % If the 'time' command writes out the results with a decimal point instead of a decimal comma, this line should be modified to: res = textscan (str, "%dm%f");
				minutes = double(res{1});
				seconds = double(res{2});
				fraction = double(res{3}); % If the 'time' command writes out the results with a decimal point instead of a decimal comma, the decimal variable should be deleted or commented out.
				total_seconds = [total_seconds; minutes * 60 + seconds + fraction/1000 ]; % If the 'time' command writes out the results with a decimal point instead of a decimal comma, this line should be modified to: total_seconds = [total_seconds; minutes * 60 + seconds ];
			end
			fclose(fid);
			q = quantile(total_seconds, [0, 0.25, 0.5, 0.75, 1]);
			stat = [stat; blocks(t), threads(t), mean(total_seconds), q];
		end
		file_save = sprintf('grid_%d_order_%d.txt', threads(1) * blocks(1), order);
		save(file_save, 'stat', '-ascii');

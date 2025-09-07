% example_animation.m - Example animation script
% Run different animation examples

fprintf('Creating Lloyd iteration animations...\n');

% % Example 1: Quick GIF with fewer stipples (faster)
% fprintf('\n1. Creating quick GIF animation (3000 stipples)...\n');
% animate_lloyd('images/example-1024px.png', 'n_stipples', 3000, 'n_iter', 50, ...
%               'format', 'gif', 'frame_delay', 0.1);

% % Example 2: High-quality MP4 
% fprintf('\n2. Creating MP4 animation (5000 stipples)...\n');
% animate_lloyd('images/example-1024px.png', 'n_stipples', 5000, 'n_iter', 50, ...
%               'format', 'mp4', 'frame_rate', 15);

% Example 3: Animation with circles instead of dots
fprintf('\n3. Creating circle-based animation...\n');
animate_lloyd('images/example-1024px.png', 'n_stipples', 10000, 'n_iter', 50, ...
              'format', 'gif', 'draw_circles', true, 'radius', 2.0);

fprintf('\nAll animations saved to stipplings/animations/ folder\n');

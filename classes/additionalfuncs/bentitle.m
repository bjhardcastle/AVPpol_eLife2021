function varargout = bentitle(titleString,textInterpreter)
if nargin < 2 || isempty( textInterpreter) ||~ischar(textInterpreter)
    textInterpreter = 'none';
end

textBoxPosition = [0 0.9 1 0.1];

fig = gcf; 
% If annotations were made previously we have to clear their text

% Return the axes that contains the annotation objects
ap = findall(fig,'Type','AnnotationPane');
for ii = 1:size(ap,2)
    for jj = 1:max(size(ap(ii),2))
        for kk = 1:max(size(ap(ii).Children))
            if isequal(ap(ii).Children(kk).Position, textBoxPosition)
                ap(ii).Children(kk).String = [];
            end
        end
    end
end

% Create new annotation with string input
h = annotation('textbox', textBoxPosition, ...
    'String', titleString , ...
    'FontWeight','bold', ...
    'EdgeColor', 'none', ...
    'Interpreter',textInterpreter,...
    'HorizontalAlignment', 'center');

if nargout
   varargout{1} = h; 
end

end
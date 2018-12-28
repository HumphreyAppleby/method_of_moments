function [base, base_prop, thedron_tri, faces, bci] = mesh_prepare(vert,thedron,face_b)
% ----- base mesh generation
[base_t, base_v, bci] = base_elm(vert,thedron,face_b);
thedron_n = size(thedron,1);
thedron_tri = zeros(thedron_n,4);
base_ccoords = zeros(thedron_n,3);

% ----- geometry properties computation
base_area = cross(vert(base_v(:,3),:)-vert(base_v(:,5),:),...
        vert(base_v(:,4),:)-vert(base_v(:,5),:));
base_area = 0.5*sqrt(base_area(:,1).^2+...
    base_area(:,2).^2+base_area(:,3).^2);

base_volume = 1/6*abs(sum((vert(thedron(:,1),:)-vert(thedron(:,4),:))...
    .*cross(vert(thedron(:,2),:)-vert(thedron(:,4),:),...
    vert(thedron(:,3),:)-vert(thedron(:,4),:)),2));

faces = [(1:size(base_v,1))' base_v(:,3:5)];
for i = 1:thedron_n
    thedron_tri(i,:) = faces(sum(ismember(faces(:,2:4),thedron(i,:)),2) == 3,1)'; 
    base_ccoords(i,:) = mean(vert(thedron(i,:),:));
end

faces = faces(:,2:4);
base_prop = [base_ccoords base_volume];
base = [base_t(:,1) base_t(:,2) base_area base_v(:,1) base_v(:,2)];
end

function [base_t, base_v, bci] = base_elm(vert,thedron,faces_b,varargin)
if (nargin == 3)
    dc = [1 0 0];
elseif (nargin == 4)
    dc = varargin{1};
else
    error('Wrong count of input arguments');
end

faceb_n = size(faces_b,1);
faces_b = [(1:faceb_n)' sort(faces_b,2)];
thedron_n = size(thedron(:,1),1);
thedron = [(1:thedron_n)' sort(thedron,2)];   
thedron_usage = zeros(thedron_n,1);
base_t = zeros(3*thedron_n,2);
base_v = zeros(3*thedron_n,5);

% ----- thedron-thedron basis functions
j = 1;
for i = 1:thedron_n
        base_left = thedron(sum(builtin('_ismemberhelper',...
        thedron((i+1):end,2:end),thedron(i,2:end)),2) == 3,1);
    if (~isempty(base_left))
        base_idx = [ones(size(base_left))*i base_left+i];
        base_t(j:(j+size(base_idx,1)-1),:) = base_idx;
        j = j+size(base_idx,1);
    end
end
roofi = j-1;

for i = 1:j-1
    thedron_usage(base_t(i,:)) = thedron_usage(base_t(i,:))+1;
    base_v(i,:) = [setxor(thedron(base_t(i,1),2:5),thedron(base_t(i,2),2:5))...
        intersect(thedron(base_t(i,1),2:5),thedron(base_t(i,2),2:5))];
    if ( dc(1)*vert(base_v(i,1),1) > dc(1)*vert(base_v(i,2),1) )
        base_v(i,:) = [base_v(i,2) base_v(i,1) base_v(i,3:5)];
        base_t(i,:) = [base_t(i,2) base_t(i,1)];
    end
end

% ----- thedron-face basis function
for i = 1:thedron_n
    if (thedron_usage(thedron(i,1),1) ~= 4)
        base_left = faces_b(sum(builtin('_ismemberhelper',faces_b(1:end,2:end),thedron(i,2:end)),2) == 3,1);
        if (~isempty(base_left))
            base_idx = [ones(size(base_left))*i base_left];
            base_t(j:(j+size(base_idx,1)-1),:) = base_idx;
            j = j+size(base_idx,1);
        end
    end
end

bci = roofi;
for i = (roofi+1):(j-1)
    roofi = roofi + 1;
    base_v(roofi,:) = [0 setxor(thedron(base_t(i,1),2:5),faces_b(base_t(i,2),2:4))...
        faces_b(base_t(i,2),2:4)];
    base_t(roofi,1) = 0;
end

base_v = base_v(1:roofi,:);
base_t = base_t(1:roofi,:);
end
function fig = mesh_result(base,base_prop,vert,faces_b,thedron,plane_coef,epoints_n)

% ----- optimization of mesh plane
object_cut = (plane_coef(1)*vert(:,1) + plane_coef(2)*vert(:,2)...
    + plane_coef(3)*vert(:,3) >= plane_coef(4)).*(1:size(vert,1))';
facesb_cut = faces_b(sum(ismember(faces_b(1:end,:),object_cut),2) <= 2 ...
    & sum(ismember(faces_b(1:end,:),object_cut),2) >= 1,:);

% ----- prepare mesh plane
vert_vec = vert(facesb_cut,1:2);
max_size = abs(max(vert_vec(:,2)));
min_size = abs(min(vert_vec(:,1)));
epoints_vec1 = linspace(-(max_size),(max_size),epoints_n);
epoints_vec2 = linspace(-(min_size),(min_size),epoints_n);
[x,y] = meshgrid(epoints_vec2,epoints_vec1);

if (plane_coef(1:3) ~= 0)
    z = (plane_coef(4) - plane_coef(1)*x - plane_coef(2)*y)/plane_coef(3);
elseif (sum(plane_coef(1:3) == 0) == 1)
    if (plane_coef(1) == 0)
        z = (plane_coef(4)-plane_coef(2)*y)/plane_coef(3);
    elseif(plane_coef(2) == 0)
        z = (plane_coef(4)-plane_coef(1)*x)/plane_coef(3);
    else
        z = x;
        x = (plane_coef(4)-plane_coef(2)*y)/plane_coef(1);
    end
elseif (sum(plane_coef(1:3) == 0) == 2)
    if (plane_coef(1) ~= 0)
        z = x;
        x = plane_coef(4)/plane_coef(1)+zeros(size(x));
    elseif (plane_coef(2) ~= 0)
        z = y;
        y = plane_coef(4)/plane_coef(2)+zeros(size(x));
    else
        z = plane_coef(4)/plane_coef(3)+zeros(size(x));
    end
end

% ----- find the centers of squares
ivert = zeros((epoints_n-1)^2,4);
X = x(1:2,1:2)-x(1:2,1);
Y = y(1:2,1:2)-y(1,1:2);
Z = z(1:2,1:2)-z(1,1);
ccoord = mean([X(1:end);Y(1:end);Z(1:end)],2)';
X = x(2:end,2:end);
Y = y(2:end,2:end);
Z = z(2:end,2:end);
ivert(:,1:3) = [X(1:end)' Y(1:end)' Z(1:end)'] - ccoord;
    
% ----- convolution of tetrahedrals with evaluation points
TR = triangulation(thedron,vert);
ivert(:,4) = TR.pointLocation(ivert(:,1:3));

mask = reshape(1:epoints_n^2,epoints_n,epoints_n);
mask = mask(1:end-1,1:end-1);
mask = mask(1:end)';
faces = [mask mask+epoints_n mask+epoints_n+1 mask+1]; 
faces = faces(~isnan(ivert(:,4)),:);

% ----- evaluation poinst computation
ivert = ivert(~isnan(ivert(:,4)),:);
ivert_n = size(ivert,1);
for k = 1:ivert_n
    evp = base(base(:,1) == ivert(k,4),:);
    evm = base(base(:,2) == ivert(k,4),:);
    evp = (vert(evp(:,4),:)-ivert(k,1:3))'*evp(:,3);
    evm = (ivert(k,1:3)-vert(evm(:,5),:))'*evm(:,3);
    ivert(k,4) = norm(evp+evm)/base_prop(ivert(k,4),4);
end

fig = newplot;
ivert(:,4) = round(63*ivert(:,4)/max(ivert(:,4))+1);
patch('Faces',faces,'Vertices',[x(:) y(:) z(:)],'FaceVertexCData',ivert(:,4),...
      'FaceColor',get(fig,'DefaultSurfaceFaceColor'), ...
      'EdgeAlpha',0.05,'Parent',fig);
end
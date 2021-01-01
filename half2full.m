function full_mesh = half2full(half_mesh_quant,data)
%HALF2FULL convert a spectral quantity from radial half mesh to radial full mesh
ns = data.ns;
full_mesh = zeros(size(half_mesh_quant));
full_mesh(:,1) = 1.5*half_mesh_quant(:,1) - 0.5*half_mesh_quant(:,2);
full_mesh(:,2:d.ns) = 0.5 * (half_mesh_quant(:,1:ns-1) + half_mesh_quant(:,2:ns));
full_mesh(:,end) = 2 * half_mesh_quant(end) - half_mesh_quant(ns-2);

end


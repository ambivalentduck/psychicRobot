function out=colorScheme(in)

out=[1 0 0;
    0 1 0;
    0 0 1;
    0 1 1;
    1 0 1;
    .96 .9 .08];

if(nargin==0)
    return
end

out=[out;rand(in-6,3)];
function q = IndexFinder(array,value)

if length(value)==1
[~, q]=min(abs(array-value));
else
    fprintf('Wrong Array Size')
    pause
end

end


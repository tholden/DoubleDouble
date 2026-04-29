function v = Assign( v, Value, varargin )
    if any( cellfun( @isempty, varargin, 'UniformOutput', true ) )
        return
    end
    v( varargin{:} ) = Value;
end

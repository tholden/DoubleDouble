function v = TimesPowerOf2( v, b )

    assert( isa( b, 'double' ) );
    v = v .* b;

end

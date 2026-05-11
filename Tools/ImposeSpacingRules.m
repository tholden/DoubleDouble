function ImposeSpacingRules()
    % ImposeSpacingRules Enforces CODING_STANDARDS.md spacing rules on all .m files

    scriptPath = fileparts( mfilename( 'fullpath' ));
    targetDir = fullfile( scriptPath, '..' );

    files = dir( fullfile( targetDir, '**', '*.m' ));

    for k = 1 : numel( files )
        filePath = fullfile( files( k ).folder, files( k ).name );

        lines = readlines( filePath );
        if isempty( lines ) || ( numel( lines ) == 1 && strlength( lines( 1 )) == 0 )
            continue;
        end

        newLines = strings( 0, 1 );

        % Phase 1: Regex replacements and collapse empty lines
        for i = 1 : numel( lines )
            line = lines( i );

            if strlength( strtrim( line )) == 0
                if i > 1 && strlength( strtrim( lines( i - 1 )) ) == 0
                    continue; % At most one empty line in a row
                end
                newLines( end + 1, 1 ) = "";
                continue;
            end

            % Split into code and comment to protect comments from regex
            idx = strfind( line, '%');
            if ~isempty( idx )
                cIdx = idx( 1 );
                codePart = extractBefore( line, cIdx );
                commentPart = extractAfter( line, cIdx - 1 );
            else
                codePart = line;
                commentPart = "";
            end

            if strlength( strtrim( codePart )) > 0
                % Rule: Open brackets followed by space
                codePart = regexprep( codePart, '( [ \[ \( \{] )( ? = [ ^\s\ ]\ )\ }] )', '$1 ' );
                % Rule: Close brackets preceded by space
                codePart = regexprep( codePart, '( [ ^\s \ [ \( \{] )( [ \ ]\ )\ }] )', '$1 $2' );

                % Rule: Colons surrounded by spaces
                codePart = regexprep( codePart, '( ? <= [ ^\s ])( : )', ' $1' );
                codePart = regexprep( codePart, '( : )( ? = [ ^\s ])', '$1 ' );

                % Rule: Commas and semicolons followed by space
                codePart = regexprep( codePart, '( [ , ; ])( ? = [ ^\s ])', '$1 ' );

                % Rule: Relational and logical operators surrounded by spaces
                codePart = regexprep( codePart, '( ? <= [ ^\s<>=~ ])( ==| ~= | <= | >= | = | < | > | && |\|\| )( ? = [ ^\s<>= ])', ' $1 ' );

                % Rule: Element-wise operators surrounded by spaces
                codePart = regexprep( codePart, '( ? <= [ ^\s ])( \ .\ *|\ .\ /|\ .\ \|\ .\ ^ )( ? = [ ^\s ])', ' $1 ' );

                % Rule: Standard operators (+, -)
                % Protect 1e-5 by avoiding negative lookbehind for e/E
                codePart = regexprep( codePart, '( ? <= [ a - zA - Z0 - 9_\ ]\ )\ }] )( ? < ![ eE ])( \+|- )( ? = [ a - zA - Z0 - 9_ \ [ \( \{] )', ' $1 ' );

                % Rule: Standard operators (*, /, \, ^, &, |)
                % Protect .* ./ etc. by avoiding negative lookbehind for .
                codePart = regexprep( codePart, '( ? <= [ a - zA - Z0 - 9_\ ]\ )\ }'''' ])( ? < !\. )( \*|\/|\\|\^|&|\| )( ? = [ a - zA - Z0 - 9_ \ [ \( \{] )', ' $1 ' );
            end

            newLines( end + 1, 1 ) = codePart + commentPart;
        end

        % Phase 2: Block spacing rules
        finalLines = strings( 0, 1 );
        for i = 1 : numel( newLines )
            finalLines( end + 1, 1 ) = newLines( i );

            % Rule: All functions, methods & properties blocks should be followed by an empty line
            match = regexp( newLines( i ), ' ^ ( \s* )( classdef | methods | properties | function ) \ b', 'tokens' );
            if ~isempty( match )
                if i < numel( newLines ) && strlength( strtrim( newLines( i + 1 )) ) > 0
                    finalLines( end + 1, 1 ) = "";
                end
            end

            % Rule: Any block indented by one or two steps should have an empty line before its end
            % (Approximated by ensuring 'end' at indent 0 or 4 spaces is preceded by empty line)
            endMatch = regexp( newLines( i ), ' ^ ( \s* )end \ b', 'tokens' );
            if ~isempty( endMatch )
                indentSpaces = endMatch{ 1 }{ 1 };
                if strlength( indentSpaces ) == 0 || strlength( indentSpaces ) == 4
                    if numel( finalLines ) > 1 && strlength( strtrim( finalLines( end - 1 )) ) > 0
                        temp = finalLines( end );
                        finalLines( end ) = "";
                        finalLines( end + 1, 1 ) = temp;
                    end
                end
            end
        end

        % Phase 3: Collapse any newly created double empty lines
        cleanLines = strings( 0, 1 );
        for i = 1 : numel( finalLines )
            if strlength( strtrim( finalLines( i )) ) == 0
                if i > 1 && strlength( strtrim( cleanLines( end )) ) == 0
                    continue;
                end
                cleanLines( end + 1, 1 ) = "";
            else
                cleanLines( end + 1, 1 ) = finalLines( i );
            end
        end

        % Write out using writelines (forces Unix line endings with '\n')
        writelines( cleanLines, filePath, 'LineEnding', ' \ n' );
    end
end

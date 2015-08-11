grammar bttest;

testfile
    : terminator? 'package' ID terminator testcase+
    ;

testcase
    : 'test' ID terminator testline+
    ;

testline
    : inputmsg terminator
    | testmsg terminator
    | enddocex terminator
    | defmat terminator
    | defvec terminator
    | defvar terminator
    | defrange terminator
    | deflist terminator
    | tcode terminator
    | tprint terminator
    | tplot terminator
    | tassert terminator
    ;

inputmsg: INPUTMSG;
testmsg: TESTMSG;
enddocex: ENDDOCEX;
defmat: 'defmat' WBRESULT? ID SEPAR '{' vecspec (',' vecspec)* '}';
defvec: 'defvec' WBRESULT? ID SEPAR vecspec;
deflist: 'deflist' WBRESULT? ID SEPAR vecspec;
defvar: 'defvar' WBRESULT? ID SEPAR expression;
defrange: 'defrange' WBRESULT? ID SEPAR NUMBER SEPAR NUMBER SEPAR NUMBER;
tplot: 'plot' ID (SEPAR ID)*;
tprint: 'print' ANNOUNCE? expression (SEPAR expression)*;
tassert: 'assert' expression SEPAR STRING;
tcode: 'code' ((ANNOUNCE? WBRESULT?) | SILENT?) (varlist '=')? expression;

vecspec: '{' expression (',' expression)* '}';

varlist
    : br='{' expression (',' expression)* '}'
    | expression
    ;

//varlist
//    : '{' ID (',' ID)* '}'
//    | ID
//    ;

expression
    :   primary
    |   tlist
    |   par1='(' expression par2=')'
    |   ID funpar='[' argumentList? ']'
    |   expression arrix='[' '[' arrindexexpr (',' arrindexexpr)* ']' ']'
    |   expression lstix='<' '<' arrindexexpr '>' '>'
    |   unop='!' expression
    |   unop='-' expression
    |   expression biop=('^'|'.^') expression
    |   expression biop=('*'|'/'|'.*'|'./'|'.') expression
    |   expression biop=('+'|'-') expression
    |   expression biop=('<=' | '>=' | '>' | '<') expression
    |   expression biop=('==' | '!=') expression
    |   expression biop='&&' expression
    |   expression biop='||' expression
    ;

arrindexexpr
    :   arrixm='-'? arrix=arrindex
    |   arrfrom=arrindex? ';' arrstep=expression? ';' ((arrendm='-'? arrend=NUMBER?) | arrende=expression)
    ;

arrindex
    :   NUMBER
    |   expression   
    ;

primary
    :   ID
    |   NUMBER
    |   STRING
    ;
    
argumentList: expression (',' expression)*;              
         
tlist
    :   '<-' expression (',' expression)* '->'
    ;
     
terminator: NEWLINE | EOF;

INPUTMSG: 'inputmsg';
TESTMSG: 'testmsg';
ENDDOCEX: 'enddocex';

NEWLINE: ('\r'? '\n')+;

WBRESULT: '@';
ANNOUNCE: '!';
SILENT: '$';
SEPAR :'#';
ID: LETTER (LETTER | [0-9] | '`')*;
fragment LETTER: [a-zA-Z];
   
STRING: '"' (~["\\])* '"';

NUMBER
    :   INT '.' INT? EXP?   // 1.35, 1.35E-9, 0.3, -4.5
    |   INT EXP            // 1e10 -3e4
    |   INT                // -3, 45
    ;

fragment INT :   [0-9]+ ;
fragment EXP :   [Ee] [+\-]? INT ; // \- since - means "range" inside [...]

WS  :   [ \t\r]+ -> skip ;
SL_COMMENT:   '//' .*? '\n' -> skip;


# Generated from java-escape by ANTLR 4.5
# encoding: utf-8
from antlr4 import *
from io import StringIO
package = globals().get("__package__", None)
ischild = len(package)>0 if package is not None else False
if ischild:
    from .bttestListener import bttestListener
else:
    from bttestListener import bttestListener
def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u0430\ud6d1\u8206\uad2d\u4417\uaef1\u8d80\uaadd\3\66")
        buf.write("\u0161\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t\7")
        buf.write("\4\b\t\b\4\t\t\t\4\n\t\n\4\13\t\13\4\f\t\f\4\r\t\r\4\16")
        buf.write("\t\16\4\17\t\17\4\20\t\20\4\21\t\21\4\22\t\22\4\23\t\23")
        buf.write("\4\24\t\24\4\25\t\25\4\26\t\26\4\27\t\27\4\30\t\30\4\31")
        buf.write("\t\31\3\2\5\2\64\n\2\3\2\3\2\3\2\3\2\6\2:\n\2\r\2\16\2")
        buf.write(";\3\3\3\3\3\3\3\3\6\3B\n\3\r\3\16\3C\3\4\3\4\3\4\3\4\3")
        buf.write("\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4")
        buf.write("\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3\4\3")
        buf.write("\4\3\4\3\4\3\4\3\4\5\4j\n\4\3\5\3\5\3\6\3\6\3\7\3\7\3")
        buf.write("\b\3\b\5\bt\n\b\3\b\3\b\3\b\3\b\3\b\3\b\7\b|\n\b\f\b\16")
        buf.write("\b\177\13\b\3\b\3\b\3\t\3\t\5\t\u0085\n\t\3\t\3\t\3\t")
        buf.write("\3\t\3\n\3\n\5\n\u008d\n\n\3\n\3\n\3\n\3\n\3\13\3\13\5")
        buf.write("\13\u0095\n\13\3\13\3\13\3\13\3\13\3\f\3\f\5\f\u009d\n")
        buf.write("\f\3\f\3\f\3\f\3\f\3\f\3\f\3\f\3\f\3\r\3\r\3\r\3\r\7\r")
        buf.write("\u00ab\n\r\f\r\16\r\u00ae\13\r\3\16\3\16\5\16\u00b2\n")
        buf.write("\16\3\16\3\16\3\16\7\16\u00b7\n\16\f\16\16\16\u00ba\13")
        buf.write("\16\3\17\3\17\3\17\3\17\3\17\3\20\3\20\5\20\u00c3\n\20")
        buf.write("\3\20\5\20\u00c6\n\20\3\20\5\20\u00c9\n\20\5\20\u00cb")
        buf.write("\n\20\3\20\3\20\3\20\5\20\u00d0\n\20\3\20\3\20\3\21\3")
        buf.write("\21\3\21\3\21\7\21\u00d8\n\21\f\21\16\21\u00db\13\21\3")
        buf.write("\21\3\21\3\22\3\22\3\22\3\22\7\22\u00e3\n\22\f\22\16\22")
        buf.write("\u00e6\13\22\3\22\3\22\3\22\5\22\u00eb\n\22\3\23\3\23")
        buf.write("\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23")
        buf.write("\3\23\5\23\u00fb\n\23\3\23\5\23\u00fe\n\23\3\23\3\23\3")
        buf.write("\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23")
        buf.write("\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23")
        buf.write("\3\23\3\23\3\23\7\23\u011b\n\23\f\23\16\23\u011e\13\23")
        buf.write("\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\3\23\7\23")
        buf.write("\u012a\n\23\f\23\16\23\u012d\13\23\3\24\5\24\u0130\n\24")
        buf.write("\3\24\3\24\5\24\u0134\n\24\3\24\3\24\5\24\u0138\n\24\3")
        buf.write("\24\3\24\5\24\u013c\n\24\3\24\5\24\u013f\n\24\3\24\5\24")
        buf.write("\u0142\n\24\5\24\u0144\n\24\3\25\3\25\5\25\u0148\n\25")
        buf.write("\3\26\3\26\3\27\3\27\3\27\7\27\u014f\n\27\f\27\16\27\u0152")
        buf.write("\13\27\3\30\3\30\3\30\3\30\7\30\u0158\n\30\f\30\16\30")
        buf.write("\u015b\13\30\3\30\3\30\3\31\3\31\3\31\2\3$\32\2\4\6\b")
        buf.write("\n\f\16\20\22\24\26\30\32\34\36 \"$&(*,.\60\2\t\3\2\31")
        buf.write("\32\3\2\33\37\4\2\30\30  \4\2\26\27!\"\3\2#$\3\2\62\64")
        buf.write("\3\3--\u0181\2\63\3\2\2\2\4=\3\2\2\2\6i\3\2\2\2\bk\3\2")
        buf.write("\2\2\nm\3\2\2\2\fo\3\2\2\2\16q\3\2\2\2\20\u0082\3\2\2")
        buf.write("\2\22\u008a\3\2\2\2\24\u0092\3\2\2\2\26\u009a\3\2\2\2")
        buf.write("\30\u00a6\3\2\2\2\32\u00af\3\2\2\2\34\u00bb\3\2\2\2\36")
        buf.write("\u00c0\3\2\2\2 \u00d3\3\2\2\2\"\u00ea\3\2\2\2$\u00fd\3")
        buf.write("\2\2\2&\u0143\3\2\2\2(\u0147\3\2\2\2*\u0149\3\2\2\2,\u014b")
        buf.write("\3\2\2\2.\u0153\3\2\2\2\60\u015e\3\2\2\2\62\64\5\60\31")
        buf.write("\2\63\62\3\2\2\2\63\64\3\2\2\2\64\65\3\2\2\2\65\66\7\3")
        buf.write("\2\2\66\67\7\62\2\2\679\5\60\31\28:\5\4\3\298\3\2\2\2")
        buf.write(":;\3\2\2\2;9\3\2\2\2;<\3\2\2\2<\3\3\2\2\2=>\7\4\2\2>?")
        buf.write("\7\62\2\2?A\5\60\31\2@B\5\6\4\2A@\3\2\2\2BC\3\2\2\2CA")
        buf.write("\3\2\2\2CD\3\2\2\2D\5\3\2\2\2EF\5\b\5\2FG\5\60\31\2Gj")
        buf.write("\3\2\2\2HI\5\n\6\2IJ\5\60\31\2Jj\3\2\2\2KL\5\f\7\2LM\5")
        buf.write("\60\31\2Mj\3\2\2\2NO\5\16\b\2OP\5\60\31\2Pj\3\2\2\2QR")
        buf.write("\5\20\t\2RS\5\60\31\2Sj\3\2\2\2TU\5\24\13\2UV\5\60\31")
        buf.write("\2Vj\3\2\2\2WX\5\26\f\2XY\5\60\31\2Yj\3\2\2\2Z[\5\22\n")
        buf.write("\2[\\\5\60\31\2\\j\3\2\2\2]^\5\36\20\2^_\5\60\31\2_j\3")
        buf.write("\2\2\2`a\5\32\16\2ab\5\60\31\2bj\3\2\2\2cd\5\30\r\2de")
        buf.write("\5\60\31\2ej\3\2\2\2fg\5\34\17\2gh\5\60\31\2hj\3\2\2\2")
        buf.write("iE\3\2\2\2iH\3\2\2\2iK\3\2\2\2iN\3\2\2\2iQ\3\2\2\2iT\3")
        buf.write("\2\2\2iW\3\2\2\2iZ\3\2\2\2i]\3\2\2\2i`\3\2\2\2ic\3\2\2")
        buf.write("\2if\3\2\2\2j\7\3\2\2\2kl\7*\2\2l\t\3\2\2\2mn\7+\2\2n")
        buf.write("\13\3\2\2\2op\7,\2\2p\r\3\2\2\2qs\7\5\2\2rt\7.\2\2sr\3")
        buf.write("\2\2\2st\3\2\2\2tu\3\2\2\2uv\7\62\2\2vw\7\61\2\2wx\7\6")
        buf.write("\2\2x}\5 \21\2yz\7\7\2\2z|\5 \21\2{y\3\2\2\2|\177\3\2")
        buf.write("\2\2}{\3\2\2\2}~\3\2\2\2~\u0080\3\2\2\2\177}\3\2\2\2\u0080")
        buf.write("\u0081\7\b\2\2\u0081\17\3\2\2\2\u0082\u0084\7\t\2\2\u0083")
        buf.write("\u0085\7.\2\2\u0084\u0083\3\2\2\2\u0084\u0085\3\2\2\2")
        buf.write("\u0085\u0086\3\2\2\2\u0086\u0087\7\62\2\2\u0087\u0088")
        buf.write("\7\61\2\2\u0088\u0089\5 \21\2\u0089\21\3\2\2\2\u008a\u008c")
        buf.write("\7\n\2\2\u008b\u008d\7.\2\2\u008c\u008b\3\2\2\2\u008c")
        buf.write("\u008d\3\2\2\2\u008d\u008e\3\2\2\2\u008e\u008f\7\62\2")
        buf.write("\2\u008f\u0090\7\61\2\2\u0090\u0091\5 \21\2\u0091\23\3")
        buf.write("\2\2\2\u0092\u0094\7\13\2\2\u0093\u0095\7.\2\2\u0094\u0093")
        buf.write("\3\2\2\2\u0094\u0095\3\2\2\2\u0095\u0096\3\2\2\2\u0096")
        buf.write("\u0097\7\62\2\2\u0097\u0098\7\61\2\2\u0098\u0099\5$\23")
        buf.write("\2\u0099\25\3\2\2\2\u009a\u009c\7\f\2\2\u009b\u009d\7")
        buf.write(".\2\2\u009c\u009b\3\2\2\2\u009c\u009d\3\2\2\2\u009d\u009e")
        buf.write("\3\2\2\2\u009e\u009f\7\62\2\2\u009f\u00a0\7\61\2\2\u00a0")
        buf.write("\u00a1\7\64\2\2\u00a1\u00a2\7\61\2\2\u00a2\u00a3\7\64")
        buf.write("\2\2\u00a3\u00a4\7\61\2\2\u00a4\u00a5\7\64\2\2\u00a5\27")
        buf.write("\3\2\2\2\u00a6\u00a7\7\r\2\2\u00a7\u00ac\7\62\2\2\u00a8")
        buf.write("\u00a9\7\61\2\2\u00a9\u00ab\7\62\2\2\u00aa\u00a8\3\2\2")
        buf.write("\2\u00ab\u00ae\3\2\2\2\u00ac\u00aa\3\2\2\2\u00ac\u00ad")
        buf.write("\3\2\2\2\u00ad\31\3\2\2\2\u00ae\u00ac\3\2\2\2\u00af\u00b1")
        buf.write("\7\16\2\2\u00b0\u00b2\7/\2\2\u00b1\u00b0\3\2\2\2\u00b1")
        buf.write("\u00b2\3\2\2\2\u00b2\u00b3\3\2\2\2\u00b3\u00b8\5$\23\2")
        buf.write("\u00b4\u00b5\7\61\2\2\u00b5\u00b7\5$\23\2\u00b6\u00b4")
        buf.write("\3\2\2\2\u00b7\u00ba\3\2\2\2\u00b8\u00b6\3\2\2\2\u00b8")
        buf.write("\u00b9\3\2\2\2\u00b9\33\3\2\2\2\u00ba\u00b8\3\2\2\2\u00bb")
        buf.write("\u00bc\7\17\2\2\u00bc\u00bd\5$\23\2\u00bd\u00be\7\61\2")
        buf.write("\2\u00be\u00bf\7\63\2\2\u00bf\35\3\2\2\2\u00c0\u00ca\7")
        buf.write("\20\2\2\u00c1\u00c3\7/\2\2\u00c2\u00c1\3\2\2\2\u00c2\u00c3")
        buf.write("\3\2\2\2\u00c3\u00c5\3\2\2\2\u00c4\u00c6\7.\2\2\u00c5")
        buf.write("\u00c4\3\2\2\2\u00c5\u00c6\3\2\2\2\u00c6\u00cb\3\2\2\2")
        buf.write("\u00c7\u00c9\7\60\2\2\u00c8\u00c7\3\2\2\2\u00c8\u00c9")
        buf.write("\3\2\2\2\u00c9\u00cb\3\2\2\2\u00ca\u00c2\3\2\2\2\u00ca")
        buf.write("\u00c8\3\2\2\2\u00cb\u00cf\3\2\2\2\u00cc\u00cd\5\"\22")
        buf.write("\2\u00cd\u00ce\7\21\2\2\u00ce\u00d0\3\2\2\2\u00cf\u00cc")
        buf.write("\3\2\2\2\u00cf\u00d0\3\2\2\2\u00d0\u00d1\3\2\2\2\u00d1")
        buf.write("\u00d2\5$\23\2\u00d2\37\3\2\2\2\u00d3\u00d4\7\6\2\2\u00d4")
        buf.write("\u00d9\5$\23\2\u00d5\u00d6\7\7\2\2\u00d6\u00d8\5$\23\2")
        buf.write("\u00d7\u00d5\3\2\2\2\u00d8\u00db\3\2\2\2\u00d9\u00d7\3")
        buf.write("\2\2\2\u00d9\u00da\3\2\2\2\u00da\u00dc\3\2\2\2\u00db\u00d9")
        buf.write("\3\2\2\2\u00dc\u00dd\7\b\2\2\u00dd!\3\2\2\2\u00de\u00df")
        buf.write("\7\6\2\2\u00df\u00e4\5$\23\2\u00e0\u00e1\7\7\2\2\u00e1")
        buf.write("\u00e3\5$\23\2\u00e2\u00e0\3\2\2\2\u00e3\u00e6\3\2\2\2")
        buf.write("\u00e4\u00e2\3\2\2\2\u00e4\u00e5\3\2\2\2\u00e5\u00e7\3")
        buf.write("\2\2\2\u00e6\u00e4\3\2\2\2\u00e7\u00e8\7\b\2\2\u00e8\u00eb")
        buf.write("\3\2\2\2\u00e9\u00eb\5$\23\2\u00ea\u00de\3\2\2\2\u00ea")
        buf.write("\u00e9\3\2\2\2\u00eb#\3\2\2\2\u00ec\u00ed\b\23\1\2\u00ed")
        buf.write("\u00ee\7/\2\2\u00ee\u00fe\5$\23\13\u00ef\u00f0\7\30\2")
        buf.write("\2\u00f0\u00fe\5$\23\n\u00f1\u00fe\5*\26\2\u00f2\u00fe")
        buf.write("\5.\30\2\u00f3\u00f4\7\22\2\2\u00f4\u00f5\5$\23\2\u00f5")
        buf.write("\u00f6\7\23\2\2\u00f6\u00fe\3\2\2\2\u00f7\u00f8\7\62\2")
        buf.write("\2\u00f8\u00fa\7\24\2\2\u00f9\u00fb\5,\27\2\u00fa\u00f9")
        buf.write("\3\2\2\2\u00fa\u00fb\3\2\2\2\u00fb\u00fc\3\2\2\2\u00fc")
        buf.write("\u00fe\7\25\2\2\u00fd\u00ec\3\2\2\2\u00fd\u00ef\3\2\2")
        buf.write("\2\u00fd\u00f1\3\2\2\2\u00fd\u00f2\3\2\2\2\u00fd\u00f3")
        buf.write("\3\2\2\2\u00fd\u00f7\3\2\2\2\u00fe\u012b\3\2\2\2\u00ff")
        buf.write("\u0100\f\t\2\2\u0100\u0101\t\2\2\2\u0101\u012a\5$\23\n")
        buf.write("\u0102\u0103\f\b\2\2\u0103\u0104\t\3\2\2\u0104\u012a\5")
        buf.write("$\23\t\u0105\u0106\f\7\2\2\u0106\u0107\t\4\2\2\u0107\u012a")
        buf.write("\5$\23\b\u0108\u0109\f\6\2\2\u0109\u010a\t\5\2\2\u010a")
        buf.write("\u012a\5$\23\7\u010b\u010c\f\5\2\2\u010c\u010d\t\6\2\2")
        buf.write("\u010d\u012a\5$\23\6\u010e\u010f\f\4\2\2\u010f\u0110\7")
        buf.write("%\2\2\u0110\u012a\5$\23\5\u0111\u0112\f\3\2\2\u0112\u0113")
        buf.write("\7&\2\2\u0113\u012a\5$\23\4\u0114\u0115\f\r\2\2\u0115")
        buf.write("\u0116\7\24\2\2\u0116\u0117\7\24\2\2\u0117\u011c\5&\24")
        buf.write("\2\u0118\u0119\7\7\2\2\u0119\u011b\5&\24\2\u011a\u0118")
        buf.write("\3\2\2\2\u011b\u011e\3\2\2\2\u011c\u011a\3\2\2\2\u011c")
        buf.write("\u011d\3\2\2\2\u011d\u011f\3\2\2\2\u011e\u011c\3\2\2\2")
        buf.write("\u011f\u0120\7\25\2\2\u0120\u0121\7\25\2\2\u0121\u012a")
        buf.write("\3\2\2\2\u0122\u0123\f\f\2\2\u0123\u0124\7\26\2\2\u0124")
        buf.write("\u0125\7\26\2\2\u0125\u0126\5&\24\2\u0126\u0127\7\27\2")
        buf.write("\2\u0127\u0128\7\27\2\2\u0128\u012a\3\2\2\2\u0129\u00ff")
        buf.write("\3\2\2\2\u0129\u0102\3\2\2\2\u0129\u0105\3\2\2\2\u0129")
        buf.write("\u0108\3\2\2\2\u0129\u010b\3\2\2\2\u0129\u010e\3\2\2\2")
        buf.write("\u0129\u0111\3\2\2\2\u0129\u0114\3\2\2\2\u0129\u0122\3")
        buf.write("\2\2\2\u012a\u012d\3\2\2\2\u012b\u0129\3\2\2\2\u012b\u012c")
        buf.write("\3\2\2\2\u012c%\3\2\2\2\u012d\u012b\3\2\2\2\u012e\u0130")
        buf.write("\7\30\2\2\u012f\u012e\3\2\2\2\u012f\u0130\3\2\2\2\u0130")
        buf.write("\u0131\3\2\2\2\u0131\u0144\5(\25\2\u0132\u0134\5(\25\2")
        buf.write("\u0133\u0132\3\2\2\2\u0133\u0134\3\2\2\2\u0134\u0135\3")
        buf.write("\2\2\2\u0135\u0137\7\'\2\2\u0136\u0138\5$\23\2\u0137\u0136")
        buf.write("\3\2\2\2\u0137\u0138\3\2\2\2\u0138\u0139\3\2\2\2\u0139")
        buf.write("\u0141\7\'\2\2\u013a\u013c\7\30\2\2\u013b\u013a\3\2\2")
        buf.write("\2\u013b\u013c\3\2\2\2\u013c\u013e\3\2\2\2\u013d\u013f")
        buf.write("\7\64\2\2\u013e\u013d\3\2\2\2\u013e\u013f\3\2\2\2\u013f")
        buf.write("\u0142\3\2\2\2\u0140\u0142\5$\23\2\u0141\u013b\3\2\2\2")
        buf.write("\u0141\u0140\3\2\2\2\u0142\u0144\3\2\2\2\u0143\u012f\3")
        buf.write("\2\2\2\u0143\u0133\3\2\2\2\u0144\'\3\2\2\2\u0145\u0148")
        buf.write("\7\64\2\2\u0146\u0148\5$\23\2\u0147\u0145\3\2\2\2\u0147")
        buf.write("\u0146\3\2\2\2\u0148)\3\2\2\2\u0149\u014a\t\7\2\2\u014a")
        buf.write("+\3\2\2\2\u014b\u0150\5$\23\2\u014c\u014d\7\7\2\2\u014d")
        buf.write("\u014f\5$\23\2\u014e\u014c\3\2\2\2\u014f\u0152\3\2\2\2")
        buf.write("\u0150\u014e\3\2\2\2\u0150\u0151\3\2\2\2\u0151-\3\2\2")
        buf.write("\2\u0152\u0150\3\2\2\2\u0153\u0154\7(\2\2\u0154\u0159")
        buf.write("\5$\23\2\u0155\u0156\7\7\2\2\u0156\u0158\5$\23\2\u0157")
        buf.write("\u0155\3\2\2\2\u0158\u015b\3\2\2\2\u0159\u0157\3\2\2\2")
        buf.write("\u0159\u015a\3\2\2\2\u015a\u015c\3\2\2\2\u015b\u0159\3")
        buf.write("\2\2\2\u015c\u015d\7)\2\2\u015d/\3\2\2\2\u015e\u015f\t")
        buf.write("\b\2\2\u015f\61\3\2\2\2&\63;Cis}\u0084\u008c\u0094\u009c")
        buf.write("\u00ac\u00b1\u00b8\u00c2\u00c5\u00c8\u00ca\u00cf\u00d9")
        buf.write("\u00e4\u00ea\u00fa\u00fd\u011c\u0129\u012b\u012f\u0133")
        buf.write("\u0137\u013b\u013e\u0141\u0143\u0147\u0150\u0159")
        return buf.getvalue()


class bttestParser ( Parser ):

    grammarFileName = "java-escape"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ u"<INVALID>", u"'package'", u"'test'", u"'defmat'", 
                     u"'{'", u"','", u"'}'", u"'defvec'", u"'deflist'", 
                     u"'defvar'", u"'defrange'", u"'plot'", u"'print'", 
                     u"'assert'", u"'code'", u"'='", u"'('", u"')'", u"'['", 
                     u"']'", u"'<'", u"'>'", u"'-'", u"'^'", u"'.^'", u"'*'", 
                     u"'/'", u"'.*'", u"'./'", u"'.'", u"'+'", u"'<='", 
                     u"'>='", u"'=='", u"'!='", u"'&&'", u"'||'", u"';'", 
                     u"'<-'", u"'->'", u"'inputmsg'", u"'testmsg'", u"'enddocex'", 
                     u"<INVALID>", u"'@'", u"'!'", u"'$'", u"'#'" ]

    symbolicNames = [ u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"<INVALID>", u"<INVALID>", u"<INVALID>", u"<INVALID>", 
                      u"INPUTMSG", u"TESTMSG", u"ENDDOCEX", u"NEWLINE", 
                      u"WBRESULT", u"ANNOUNCE", u"SILENT", u"SEPAR", u"ID", 
                      u"STRING", u"NUMBER", u"WS", u"SL_COMMENT" ]

    RULE_testfile = 0
    RULE_testcase = 1
    RULE_testline = 2
    RULE_inputmsg = 3
    RULE_testmsg = 4
    RULE_enddocex = 5
    RULE_defmat = 6
    RULE_defvec = 7
    RULE_deflist = 8
    RULE_defvar = 9
    RULE_defrange = 10
    RULE_tplot = 11
    RULE_tprint = 12
    RULE_tassert = 13
    RULE_tcode = 14
    RULE_vecspec = 15
    RULE_varlist = 16
    RULE_expression = 17
    RULE_arrindexexpr = 18
    RULE_arrindex = 19
    RULE_primary = 20
    RULE_argumentList = 21
    RULE_tlist = 22
    RULE_terminator = 23

    ruleNames =  [ "testfile", "testcase", "testline", "inputmsg", "testmsg", 
                   "enddocex", "defmat", "defvec", "deflist", "defvar", 
                   "defrange", "tplot", "tprint", "tassert", "tcode", "vecspec", 
                   "varlist", "expression", "arrindexexpr", "arrindex", 
                   "primary", "argumentList", "tlist", "terminator" ]

    EOF = Token.EOF
    T__0=1
    T__1=2
    T__2=3
    T__3=4
    T__4=5
    T__5=6
    T__6=7
    T__7=8
    T__8=9
    T__9=10
    T__10=11
    T__11=12
    T__12=13
    T__13=14
    T__14=15
    T__15=16
    T__16=17
    T__17=18
    T__18=19
    T__19=20
    T__20=21
    T__21=22
    T__22=23
    T__23=24
    T__24=25
    T__25=26
    T__26=27
    T__27=28
    T__28=29
    T__29=30
    T__30=31
    T__31=32
    T__32=33
    T__33=34
    T__34=35
    T__35=36
    T__36=37
    T__37=38
    T__38=39
    INPUTMSG=40
    TESTMSG=41
    ENDDOCEX=42
    NEWLINE=43
    WBRESULT=44
    ANNOUNCE=45
    SILENT=46
    SEPAR=47
    ID=48
    STRING=49
    NUMBER=50
    WS=51
    SL_COMMENT=52

    def __init__(self, input:TokenStream):
        super().__init__(input)
        self.checkVersion("4.5")
        self._interp = ParserATNSimulator(self, self.atn, self.decisionsToDFA, self.sharedContextCache)
        self._predicates = None



    class TestfileContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def terminator(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.TerminatorContext)
            else:
                return self.getTypedRuleContext(bttestParser.TerminatorContext,i)


        def testcase(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.TestcaseContext)
            else:
                return self.getTypedRuleContext(bttestParser.TestcaseContext,i)


        def getRuleIndex(self):
            return bttestParser.RULE_testfile

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTestfile(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTestfile(self)




    def testfile(self):

        localctx = bttestParser.TestfileContext(self, self._ctx, self.state)
        self.enterRule(localctx, 0, self.RULE_testfile)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 49
            _la = self._input.LA(1)
            if _la==bttestParser.EOF or _la==bttestParser.NEWLINE:
                self.state = 48
                self.terminator()


            self.state = 51
            self.match(bttestParser.T__0)
            self.state = 52
            self.match(bttestParser.ID)
            self.state = 53
            self.terminator()
            self.state = 55 
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while True:
                self.state = 54
                self.testcase()
                self.state = 57 
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if not (_la==bttestParser.T__1):
                    break

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TestcaseContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def terminator(self):
            return self.getTypedRuleContext(bttestParser.TerminatorContext,0)


        def testline(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.TestlineContext)
            else:
                return self.getTypedRuleContext(bttestParser.TestlineContext,i)


        def getRuleIndex(self):
            return bttestParser.RULE_testcase

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTestcase(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTestcase(self)




    def testcase(self):

        localctx = bttestParser.TestcaseContext(self, self._ctx, self.state)
        self.enterRule(localctx, 2, self.RULE_testcase)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 59
            self.match(bttestParser.T__1)
            self.state = 60
            self.match(bttestParser.ID)
            self.state = 61
            self.terminator()
            self.state = 63 
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while True:
                self.state = 62
                self.testline()
                self.state = 65 
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if not ((((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << bttestParser.T__2) | (1 << bttestParser.T__6) | (1 << bttestParser.T__7) | (1 << bttestParser.T__8) | (1 << bttestParser.T__9) | (1 << bttestParser.T__10) | (1 << bttestParser.T__11) | (1 << bttestParser.T__12) | (1 << bttestParser.T__13) | (1 << bttestParser.INPUTMSG) | (1 << bttestParser.TESTMSG) | (1 << bttestParser.ENDDOCEX))) != 0)):
                    break

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TestlineContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def inputmsg(self):
            return self.getTypedRuleContext(bttestParser.InputmsgContext,0)


        def terminator(self):
            return self.getTypedRuleContext(bttestParser.TerminatorContext,0)


        def testmsg(self):
            return self.getTypedRuleContext(bttestParser.TestmsgContext,0)


        def enddocex(self):
            return self.getTypedRuleContext(bttestParser.EnddocexContext,0)


        def defmat(self):
            return self.getTypedRuleContext(bttestParser.DefmatContext,0)


        def defvec(self):
            return self.getTypedRuleContext(bttestParser.DefvecContext,0)


        def defvar(self):
            return self.getTypedRuleContext(bttestParser.DefvarContext,0)


        def defrange(self):
            return self.getTypedRuleContext(bttestParser.DefrangeContext,0)


        def deflist(self):
            return self.getTypedRuleContext(bttestParser.DeflistContext,0)


        def tcode(self):
            return self.getTypedRuleContext(bttestParser.TcodeContext,0)


        def tprint(self):
            return self.getTypedRuleContext(bttestParser.TprintContext,0)


        def tplot(self):
            return self.getTypedRuleContext(bttestParser.TplotContext,0)


        def tassert(self):
            return self.getTypedRuleContext(bttestParser.TassertContext,0)


        def getRuleIndex(self):
            return bttestParser.RULE_testline

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTestline(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTestline(self)




    def testline(self):

        localctx = bttestParser.TestlineContext(self, self._ctx, self.state)
        self.enterRule(localctx, 4, self.RULE_testline)
        try:
            self.state = 103
            token = self._input.LA(1)
            if token in [bttestParser.INPUTMSG]:
                self.enterOuterAlt(localctx, 1)
                self.state = 67
                self.inputmsg()
                self.state = 68
                self.terminator()

            elif token in [bttestParser.TESTMSG]:
                self.enterOuterAlt(localctx, 2)
                self.state = 70
                self.testmsg()
                self.state = 71
                self.terminator()

            elif token in [bttestParser.ENDDOCEX]:
                self.enterOuterAlt(localctx, 3)
                self.state = 73
                self.enddocex()
                self.state = 74
                self.terminator()

            elif token in [bttestParser.T__2]:
                self.enterOuterAlt(localctx, 4)
                self.state = 76
                self.defmat()
                self.state = 77
                self.terminator()

            elif token in [bttestParser.T__6]:
                self.enterOuterAlt(localctx, 5)
                self.state = 79
                self.defvec()
                self.state = 80
                self.terminator()

            elif token in [bttestParser.T__8]:
                self.enterOuterAlt(localctx, 6)
                self.state = 82
                self.defvar()
                self.state = 83
                self.terminator()

            elif token in [bttestParser.T__9]:
                self.enterOuterAlt(localctx, 7)
                self.state = 85
                self.defrange()
                self.state = 86
                self.terminator()

            elif token in [bttestParser.T__7]:
                self.enterOuterAlt(localctx, 8)
                self.state = 88
                self.deflist()
                self.state = 89
                self.terminator()

            elif token in [bttestParser.T__13]:
                self.enterOuterAlt(localctx, 9)
                self.state = 91
                self.tcode()
                self.state = 92
                self.terminator()

            elif token in [bttestParser.T__11]:
                self.enterOuterAlt(localctx, 10)
                self.state = 94
                self.tprint()
                self.state = 95
                self.terminator()

            elif token in [bttestParser.T__10]:
                self.enterOuterAlt(localctx, 11)
                self.state = 97
                self.tplot()
                self.state = 98
                self.terminator()

            elif token in [bttestParser.T__12]:
                self.enterOuterAlt(localctx, 12)
                self.state = 100
                self.tassert()
                self.state = 101
                self.terminator()

            else:
                raise NoViableAltException(self)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class InputmsgContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def INPUTMSG(self):
            return self.getToken(bttestParser.INPUTMSG, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_inputmsg

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterInputmsg(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitInputmsg(self)




    def inputmsg(self):

        localctx = bttestParser.InputmsgContext(self, self._ctx, self.state)
        self.enterRule(localctx, 6, self.RULE_inputmsg)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 105
            self.match(bttestParser.INPUTMSG)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TestmsgContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def TESTMSG(self):
            return self.getToken(bttestParser.TESTMSG, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_testmsg

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTestmsg(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTestmsg(self)




    def testmsg(self):

        localctx = bttestParser.TestmsgContext(self, self._ctx, self.state)
        self.enterRule(localctx, 8, self.RULE_testmsg)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 107
            self.match(bttestParser.TESTMSG)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class EnddocexContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ENDDOCEX(self):
            return self.getToken(bttestParser.ENDDOCEX, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_enddocex

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterEnddocex(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitEnddocex(self)




    def enddocex(self):

        localctx = bttestParser.EnddocexContext(self, self._ctx, self.state)
        self.enterRule(localctx, 10, self.RULE_enddocex)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 109
            self.match(bttestParser.ENDDOCEX)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class DefmatContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def SEPAR(self):
            return self.getToken(bttestParser.SEPAR, 0)

        def vecspec(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.VecspecContext)
            else:
                return self.getTypedRuleContext(bttestParser.VecspecContext,i)


        def WBRESULT(self):
            return self.getToken(bttestParser.WBRESULT, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_defmat

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterDefmat(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitDefmat(self)




    def defmat(self):

        localctx = bttestParser.DefmatContext(self, self._ctx, self.state)
        self.enterRule(localctx, 12, self.RULE_defmat)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 111
            self.match(bttestParser.T__2)
            self.state = 113
            _la = self._input.LA(1)
            if _la==bttestParser.WBRESULT:
                self.state = 112
                self.match(bttestParser.WBRESULT)


            self.state = 115
            self.match(bttestParser.ID)
            self.state = 116
            self.match(bttestParser.SEPAR)
            self.state = 117
            self.match(bttestParser.T__3)
            self.state = 118
            self.vecspec()
            self.state = 123
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==bttestParser.T__4:
                self.state = 119
                self.match(bttestParser.T__4)
                self.state = 120
                self.vecspec()
                self.state = 125
                self._errHandler.sync(self)
                _la = self._input.LA(1)

            self.state = 126
            self.match(bttestParser.T__5)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class DefvecContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def SEPAR(self):
            return self.getToken(bttestParser.SEPAR, 0)

        def vecspec(self):
            return self.getTypedRuleContext(bttestParser.VecspecContext,0)


        def WBRESULT(self):
            return self.getToken(bttestParser.WBRESULT, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_defvec

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterDefvec(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitDefvec(self)




    def defvec(self):

        localctx = bttestParser.DefvecContext(self, self._ctx, self.state)
        self.enterRule(localctx, 14, self.RULE_defvec)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 128
            self.match(bttestParser.T__6)
            self.state = 130
            _la = self._input.LA(1)
            if _la==bttestParser.WBRESULT:
                self.state = 129
                self.match(bttestParser.WBRESULT)


            self.state = 132
            self.match(bttestParser.ID)
            self.state = 133
            self.match(bttestParser.SEPAR)
            self.state = 134
            self.vecspec()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class DeflistContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def SEPAR(self):
            return self.getToken(bttestParser.SEPAR, 0)

        def vecspec(self):
            return self.getTypedRuleContext(bttestParser.VecspecContext,0)


        def WBRESULT(self):
            return self.getToken(bttestParser.WBRESULT, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_deflist

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterDeflist(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitDeflist(self)




    def deflist(self):

        localctx = bttestParser.DeflistContext(self, self._ctx, self.state)
        self.enterRule(localctx, 16, self.RULE_deflist)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 136
            self.match(bttestParser.T__7)
            self.state = 138
            _la = self._input.LA(1)
            if _la==bttestParser.WBRESULT:
                self.state = 137
                self.match(bttestParser.WBRESULT)


            self.state = 140
            self.match(bttestParser.ID)
            self.state = 141
            self.match(bttestParser.SEPAR)
            self.state = 142
            self.vecspec()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class DefvarContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def SEPAR(self):
            return self.getToken(bttestParser.SEPAR, 0)

        def expression(self):
            return self.getTypedRuleContext(bttestParser.ExpressionContext,0)


        def WBRESULT(self):
            return self.getToken(bttestParser.WBRESULT, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_defvar

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterDefvar(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitDefvar(self)




    def defvar(self):

        localctx = bttestParser.DefvarContext(self, self._ctx, self.state)
        self.enterRule(localctx, 18, self.RULE_defvar)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 144
            self.match(bttestParser.T__8)
            self.state = 146
            _la = self._input.LA(1)
            if _la==bttestParser.WBRESULT:
                self.state = 145
                self.match(bttestParser.WBRESULT)


            self.state = 148
            self.match(bttestParser.ID)
            self.state = 149
            self.match(bttestParser.SEPAR)
            self.state = 150
            self.expression(0)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class DefrangeContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def SEPAR(self, i:int=None):
            if i is None:
                return self.getTokens(bttestParser.SEPAR)
            else:
                return self.getToken(bttestParser.SEPAR, i)

        def NUMBER(self, i:int=None):
            if i is None:
                return self.getTokens(bttestParser.NUMBER)
            else:
                return self.getToken(bttestParser.NUMBER, i)

        def WBRESULT(self):
            return self.getToken(bttestParser.WBRESULT, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_defrange

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterDefrange(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitDefrange(self)




    def defrange(self):

        localctx = bttestParser.DefrangeContext(self, self._ctx, self.state)
        self.enterRule(localctx, 20, self.RULE_defrange)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 152
            self.match(bttestParser.T__9)
            self.state = 154
            _la = self._input.LA(1)
            if _la==bttestParser.WBRESULT:
                self.state = 153
                self.match(bttestParser.WBRESULT)


            self.state = 156
            self.match(bttestParser.ID)
            self.state = 157
            self.match(bttestParser.SEPAR)
            self.state = 158
            self.match(bttestParser.NUMBER)
            self.state = 159
            self.match(bttestParser.SEPAR)
            self.state = 160
            self.match(bttestParser.NUMBER)
            self.state = 161
            self.match(bttestParser.SEPAR)
            self.state = 162
            self.match(bttestParser.NUMBER)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TplotContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self, i:int=None):
            if i is None:
                return self.getTokens(bttestParser.ID)
            else:
                return self.getToken(bttestParser.ID, i)

        def SEPAR(self, i:int=None):
            if i is None:
                return self.getTokens(bttestParser.SEPAR)
            else:
                return self.getToken(bttestParser.SEPAR, i)

        def getRuleIndex(self):
            return bttestParser.RULE_tplot

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTplot(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTplot(self)




    def tplot(self):

        localctx = bttestParser.TplotContext(self, self._ctx, self.state)
        self.enterRule(localctx, 22, self.RULE_tplot)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 164
            self.match(bttestParser.T__10)
            self.state = 165
            self.match(bttestParser.ID)
            self.state = 170
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==bttestParser.SEPAR:
                self.state = 166
                self.match(bttestParser.SEPAR)
                self.state = 167
                self.match(bttestParser.ID)
                self.state = 172
                self._errHandler.sync(self)
                _la = self._input.LA(1)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TprintContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def expression(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.ExpressionContext)
            else:
                return self.getTypedRuleContext(bttestParser.ExpressionContext,i)


        def ANNOUNCE(self):
            return self.getToken(bttestParser.ANNOUNCE, 0)

        def SEPAR(self, i:int=None):
            if i is None:
                return self.getTokens(bttestParser.SEPAR)
            else:
                return self.getToken(bttestParser.SEPAR, i)

        def getRuleIndex(self):
            return bttestParser.RULE_tprint

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTprint(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTprint(self)




    def tprint(self):

        localctx = bttestParser.TprintContext(self, self._ctx, self.state)
        self.enterRule(localctx, 24, self.RULE_tprint)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 173
            self.match(bttestParser.T__11)
            self.state = 175
            la_ = self._interp.adaptivePredict(self._input,11,self._ctx)
            if la_ == 1:
                self.state = 174
                self.match(bttestParser.ANNOUNCE)


            self.state = 177
            self.expression(0)
            self.state = 182
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==bttestParser.SEPAR:
                self.state = 178
                self.match(bttestParser.SEPAR)
                self.state = 179
                self.expression(0)
                self.state = 184
                self._errHandler.sync(self)
                _la = self._input.LA(1)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TassertContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def expression(self):
            return self.getTypedRuleContext(bttestParser.ExpressionContext,0)


        def SEPAR(self):
            return self.getToken(bttestParser.SEPAR, 0)

        def STRING(self):
            return self.getToken(bttestParser.STRING, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_tassert

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTassert(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTassert(self)




    def tassert(self):

        localctx = bttestParser.TassertContext(self, self._ctx, self.state)
        self.enterRule(localctx, 26, self.RULE_tassert)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 185
            self.match(bttestParser.T__12)
            self.state = 186
            self.expression(0)
            self.state = 187
            self.match(bttestParser.SEPAR)
            self.state = 188
            self.match(bttestParser.STRING)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TcodeContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def expression(self):
            return self.getTypedRuleContext(bttestParser.ExpressionContext,0)


        def varlist(self):
            return self.getTypedRuleContext(bttestParser.VarlistContext,0)


        def SILENT(self):
            return self.getToken(bttestParser.SILENT, 0)

        def ANNOUNCE(self):
            return self.getToken(bttestParser.ANNOUNCE, 0)

        def WBRESULT(self):
            return self.getToken(bttestParser.WBRESULT, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_tcode

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTcode(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTcode(self)




    def tcode(self):

        localctx = bttestParser.TcodeContext(self, self._ctx, self.state)
        self.enterRule(localctx, 28, self.RULE_tcode)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 190
            self.match(bttestParser.T__13)
            self.state = 200
            la_ = self._interp.adaptivePredict(self._input,16,self._ctx)
            if la_ == 1:
                self.state = 192
                la_ = self._interp.adaptivePredict(self._input,13,self._ctx)
                if la_ == 1:
                    self.state = 191
                    self.match(bttestParser.ANNOUNCE)


                self.state = 195
                _la = self._input.LA(1)
                if _la==bttestParser.WBRESULT:
                    self.state = 194
                    self.match(bttestParser.WBRESULT)


                pass

            elif la_ == 2:
                self.state = 198
                _la = self._input.LA(1)
                if _la==bttestParser.SILENT:
                    self.state = 197
                    self.match(bttestParser.SILENT)


                pass


            self.state = 205
            la_ = self._interp.adaptivePredict(self._input,17,self._ctx)
            if la_ == 1:
                self.state = 202
                self.varlist()
                self.state = 203
                self.match(bttestParser.T__14)


            self.state = 207
            self.expression(0)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class VecspecContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def expression(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.ExpressionContext)
            else:
                return self.getTypedRuleContext(bttestParser.ExpressionContext,i)


        def getRuleIndex(self):
            return bttestParser.RULE_vecspec

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterVecspec(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitVecspec(self)




    def vecspec(self):

        localctx = bttestParser.VecspecContext(self, self._ctx, self.state)
        self.enterRule(localctx, 30, self.RULE_vecspec)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 209
            self.match(bttestParser.T__3)
            self.state = 210
            self.expression(0)
            self.state = 215
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==bttestParser.T__4:
                self.state = 211
                self.match(bttestParser.T__4)
                self.state = 212
                self.expression(0)
                self.state = 217
                self._errHandler.sync(self)
                _la = self._input.LA(1)

            self.state = 218
            self.match(bttestParser.T__5)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class VarlistContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser
            self.br = None # Token

        def expression(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.ExpressionContext)
            else:
                return self.getTypedRuleContext(bttestParser.ExpressionContext,i)


        def getRuleIndex(self):
            return bttestParser.RULE_varlist

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterVarlist(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitVarlist(self)




    def varlist(self):

        localctx = bttestParser.VarlistContext(self, self._ctx, self.state)
        self.enterRule(localctx, 32, self.RULE_varlist)
        self._la = 0 # Token type
        try:
            self.state = 232
            token = self._input.LA(1)
            if token in [bttestParser.T__3]:
                self.enterOuterAlt(localctx, 1)
                self.state = 220
                localctx.br = self.match(bttestParser.T__3)
                self.state = 221
                self.expression(0)
                self.state = 226
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==bttestParser.T__4:
                    self.state = 222
                    self.match(bttestParser.T__4)
                    self.state = 223
                    self.expression(0)
                    self.state = 228
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                self.state = 229
                self.match(bttestParser.T__5)

            elif token in [bttestParser.T__15, bttestParser.T__21, bttestParser.T__37, bttestParser.ANNOUNCE, bttestParser.ID, bttestParser.STRING, bttestParser.NUMBER]:
                self.enterOuterAlt(localctx, 2)
                self.state = 231
                self.expression(0)

            else:
                raise NoViableAltException(self)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class ExpressionContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser
            self.unop = None # Token
            self.par1 = None # Token
            self.par2 = None # Token
            self.funpar = None # Token
            self.biop = None # Token
            self.arrix = None # Token
            self.lstix = None # Token

        def expression(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.ExpressionContext)
            else:
                return self.getTypedRuleContext(bttestParser.ExpressionContext,i)


        def primary(self):
            return self.getTypedRuleContext(bttestParser.PrimaryContext,0)


        def tlist(self):
            return self.getTypedRuleContext(bttestParser.TlistContext,0)


        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def argumentList(self):
            return self.getTypedRuleContext(bttestParser.ArgumentListContext,0)


        def arrindexexpr(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.ArrindexexprContext)
            else:
                return self.getTypedRuleContext(bttestParser.ArrindexexprContext,i)


        def getRuleIndex(self):
            return bttestParser.RULE_expression

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterExpression(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitExpression(self)



    def expression(self, _p:int=0):
        _parentctx = self._ctx
        _parentState = self.state
        localctx = bttestParser.ExpressionContext(self, self._ctx, _parentState)
        _prevctx = localctx
        _startState = 34
        self.enterRecursionRule(localctx, 34, self.RULE_expression, _p)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 251
            la_ = self._interp.adaptivePredict(self._input,22,self._ctx)
            if la_ == 1:
                self.state = 235
                localctx.unop = self.match(bttestParser.ANNOUNCE)
                self.state = 236
                self.expression(9)
                pass

            elif la_ == 2:
                self.state = 237
                localctx.unop = self.match(bttestParser.T__21)
                self.state = 238
                self.expression(8)
                pass

            elif la_ == 3:
                self.state = 239
                self.primary()
                pass

            elif la_ == 4:
                self.state = 240
                self.tlist()
                pass

            elif la_ == 5:
                self.state = 241
                localctx.par1 = self.match(bttestParser.T__15)
                self.state = 242
                self.expression(0)
                self.state = 243
                localctx.par2 = self.match(bttestParser.T__16)
                pass

            elif la_ == 6:
                self.state = 245
                self.match(bttestParser.ID)
                self.state = 246
                localctx.funpar = self.match(bttestParser.T__17)
                self.state = 248
                _la = self._input.LA(1)
                if (((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << bttestParser.T__15) | (1 << bttestParser.T__21) | (1 << bttestParser.T__37) | (1 << bttestParser.ANNOUNCE) | (1 << bttestParser.ID) | (1 << bttestParser.STRING) | (1 << bttestParser.NUMBER))) != 0):
                    self.state = 247
                    self.argumentList()


                self.state = 250
                self.match(bttestParser.T__18)
                pass


            self._ctx.stop = self._input.LT(-1)
            self.state = 297
            self._errHandler.sync(self)
            _alt = self._interp.adaptivePredict(self._input,25,self._ctx)
            while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                if _alt==1:
                    if self._parseListeners is not None:
                        self.triggerExitRuleEvent()
                    _prevctx = localctx
                    self.state = 295
                    la_ = self._interp.adaptivePredict(self._input,24,self._ctx)
                    if la_ == 1:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 253
                        if not self.precpred(self._ctx, 7):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 7)")
                        self.state = 254
                        localctx.biop = self._input.LT(1)
                        _la = self._input.LA(1)
                        if not(_la==bttestParser.T__22 or _la==bttestParser.T__23):
                            localctx.biop = self._errHandler.recoverInline(self)
                        else:
                            self.consume()
                        self.state = 255
                        self.expression(8)
                        pass

                    elif la_ == 2:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 256
                        if not self.precpred(self._ctx, 6):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 6)")
                        self.state = 257
                        localctx.biop = self._input.LT(1)
                        _la = self._input.LA(1)
                        if not((((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << bttestParser.T__24) | (1 << bttestParser.T__25) | (1 << bttestParser.T__26) | (1 << bttestParser.T__27) | (1 << bttestParser.T__28))) != 0)):
                            localctx.biop = self._errHandler.recoverInline(self)
                        else:
                            self.consume()
                        self.state = 258
                        self.expression(7)
                        pass

                    elif la_ == 3:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 259
                        if not self.precpred(self._ctx, 5):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 5)")
                        self.state = 260
                        localctx.biop = self._input.LT(1)
                        _la = self._input.LA(1)
                        if not(_la==bttestParser.T__21 or _la==bttestParser.T__29):
                            localctx.biop = self._errHandler.recoverInline(self)
                        else:
                            self.consume()
                        self.state = 261
                        self.expression(6)
                        pass

                    elif la_ == 4:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 262
                        if not self.precpred(self._ctx, 4):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 4)")
                        self.state = 263
                        localctx.biop = self._input.LT(1)
                        _la = self._input.LA(1)
                        if not((((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << bttestParser.T__19) | (1 << bttestParser.T__20) | (1 << bttestParser.T__30) | (1 << bttestParser.T__31))) != 0)):
                            localctx.biop = self._errHandler.recoverInline(self)
                        else:
                            self.consume()
                        self.state = 264
                        self.expression(5)
                        pass

                    elif la_ == 5:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 265
                        if not self.precpred(self._ctx, 3):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 3)")
                        self.state = 266
                        localctx.biop = self._input.LT(1)
                        _la = self._input.LA(1)
                        if not(_la==bttestParser.T__32 or _la==bttestParser.T__33):
                            localctx.biop = self._errHandler.recoverInline(self)
                        else:
                            self.consume()
                        self.state = 267
                        self.expression(4)
                        pass

                    elif la_ == 6:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 268
                        if not self.precpred(self._ctx, 2):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 2)")
                        self.state = 269
                        localctx.biop = self.match(bttestParser.T__34)
                        self.state = 270
                        self.expression(3)
                        pass

                    elif la_ == 7:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 271
                        if not self.precpred(self._ctx, 1):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 1)")
                        self.state = 272
                        localctx.biop = self.match(bttestParser.T__35)
                        self.state = 273
                        self.expression(2)
                        pass

                    elif la_ == 8:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 274
                        if not self.precpred(self._ctx, 11):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 11)")
                        self.state = 275
                        localctx.arrix = self.match(bttestParser.T__17)
                        self.state = 276
                        self.match(bttestParser.T__17)
                        self.state = 277
                        self.arrindexexpr()
                        self.state = 282
                        self._errHandler.sync(self)
                        _la = self._input.LA(1)
                        while _la==bttestParser.T__4:
                            self.state = 278
                            self.match(bttestParser.T__4)
                            self.state = 279
                            self.arrindexexpr()
                            self.state = 284
                            self._errHandler.sync(self)
                            _la = self._input.LA(1)

                        self.state = 285
                        self.match(bttestParser.T__18)
                        self.state = 286
                        self.match(bttestParser.T__18)
                        pass

                    elif la_ == 9:
                        localctx = bttestParser.ExpressionContext(self, _parentctx, _parentState)
                        self.pushNewRecursionContext(localctx, _startState, self.RULE_expression)
                        self.state = 288
                        if not self.precpred(self._ctx, 10):
                            from antlr4.error.Errors import FailedPredicateException
                            raise FailedPredicateException(self, "self.precpred(self._ctx, 10)")
                        self.state = 289
                        localctx.lstix = self.match(bttestParser.T__19)
                        self.state = 290
                        self.match(bttestParser.T__19)
                        self.state = 291
                        self.arrindexexpr()
                        self.state = 292
                        self.match(bttestParser.T__20)
                        self.state = 293
                        self.match(bttestParser.T__20)
                        pass

             
                self.state = 299
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,25,self._ctx)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.unrollRecursionContexts(_parentctx)
        return localctx

    class ArrindexexprContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser
            self.arrixm = None # Token
            self.arrix = None # ArrindexContext
            self.arrfrom = None # ArrindexContext
            self.arrstep = None # ExpressionContext
            self.arrendm = None # Token
            self.arrend = None # Token
            self.arrende = None # ExpressionContext

        def arrindex(self):
            return self.getTypedRuleContext(bttestParser.ArrindexContext,0)


        def expression(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.ExpressionContext)
            else:
                return self.getTypedRuleContext(bttestParser.ExpressionContext,i)


        def NUMBER(self):
            return self.getToken(bttestParser.NUMBER, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_arrindexexpr

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterArrindexexpr(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitArrindexexpr(self)




    def arrindexexpr(self):

        localctx = bttestParser.ArrindexexprContext(self, self._ctx, self.state)
        self.enterRule(localctx, 36, self.RULE_arrindexexpr)
        self._la = 0 # Token type
        try:
            self.state = 321
            la_ = self._interp.adaptivePredict(self._input,32,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 301
                la_ = self._interp.adaptivePredict(self._input,26,self._ctx)
                if la_ == 1:
                    self.state = 300
                    localctx.arrixm = self.match(bttestParser.T__21)


                self.state = 303
                localctx.arrix = self.arrindex()
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 305
                _la = self._input.LA(1)
                if (((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << bttestParser.T__15) | (1 << bttestParser.T__21) | (1 << bttestParser.T__37) | (1 << bttestParser.ANNOUNCE) | (1 << bttestParser.ID) | (1 << bttestParser.STRING) | (1 << bttestParser.NUMBER))) != 0):
                    self.state = 304
                    localctx.arrfrom = self.arrindex()


                self.state = 307
                self.match(bttestParser.T__36)
                self.state = 309
                _la = self._input.LA(1)
                if (((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << bttestParser.T__15) | (1 << bttestParser.T__21) | (1 << bttestParser.T__37) | (1 << bttestParser.ANNOUNCE) | (1 << bttestParser.ID) | (1 << bttestParser.STRING) | (1 << bttestParser.NUMBER))) != 0):
                    self.state = 308
                    localctx.arrstep = self.expression(0)


                self.state = 311
                self.match(bttestParser.T__36)
                self.state = 319
                la_ = self._interp.adaptivePredict(self._input,31,self._ctx)
                if la_ == 1:
                    self.state = 313
                    _la = self._input.LA(1)
                    if _la==bttestParser.T__21:
                        self.state = 312
                        localctx.arrendm = self.match(bttestParser.T__21)


                    self.state = 316
                    _la = self._input.LA(1)
                    if _la==bttestParser.NUMBER:
                        self.state = 315
                        localctx.arrend = self.match(bttestParser.NUMBER)


                    pass

                elif la_ == 2:
                    self.state = 318
                    localctx.arrende = self.expression(0)
                    pass


                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class ArrindexContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def NUMBER(self):
            return self.getToken(bttestParser.NUMBER, 0)

        def expression(self):
            return self.getTypedRuleContext(bttestParser.ExpressionContext,0)


        def getRuleIndex(self):
            return bttestParser.RULE_arrindex

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterArrindex(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitArrindex(self)




    def arrindex(self):

        localctx = bttestParser.ArrindexContext(self, self._ctx, self.state)
        self.enterRule(localctx, 38, self.RULE_arrindex)
        try:
            self.state = 325
            la_ = self._interp.adaptivePredict(self._input,33,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 323
                self.match(bttestParser.NUMBER)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 324
                self.expression(0)
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class PrimaryContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def ID(self):
            return self.getToken(bttestParser.ID, 0)

        def NUMBER(self):
            return self.getToken(bttestParser.NUMBER, 0)

        def STRING(self):
            return self.getToken(bttestParser.STRING, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_primary

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterPrimary(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitPrimary(self)




    def primary(self):

        localctx = bttestParser.PrimaryContext(self, self._ctx, self.state)
        self.enterRule(localctx, 40, self.RULE_primary)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 327
            _la = self._input.LA(1)
            if not((((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << bttestParser.ID) | (1 << bttestParser.STRING) | (1 << bttestParser.NUMBER))) != 0)):
                self._errHandler.recoverInline(self)
            else:
                self.consume()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class ArgumentListContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def expression(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.ExpressionContext)
            else:
                return self.getTypedRuleContext(bttestParser.ExpressionContext,i)


        def getRuleIndex(self):
            return bttestParser.RULE_argumentList

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterArgumentList(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitArgumentList(self)




    def argumentList(self):

        localctx = bttestParser.ArgumentListContext(self, self._ctx, self.state)
        self.enterRule(localctx, 42, self.RULE_argumentList)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 329
            self.expression(0)
            self.state = 334
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==bttestParser.T__4:
                self.state = 330
                self.match(bttestParser.T__4)
                self.state = 331
                self.expression(0)
                self.state = 336
                self._errHandler.sync(self)
                _la = self._input.LA(1)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TlistContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def expression(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(bttestParser.ExpressionContext)
            else:
                return self.getTypedRuleContext(bttestParser.ExpressionContext,i)


        def getRuleIndex(self):
            return bttestParser.RULE_tlist

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTlist(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTlist(self)




    def tlist(self):

        localctx = bttestParser.TlistContext(self, self._ctx, self.state)
        self.enterRule(localctx, 44, self.RULE_tlist)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 337
            self.match(bttestParser.T__37)
            self.state = 338
            self.expression(0)
            self.state = 343
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==bttestParser.T__4:
                self.state = 339
                self.match(bttestParser.T__4)
                self.state = 340
                self.expression(0)
                self.state = 345
                self._errHandler.sync(self)
                _la = self._input.LA(1)

            self.state = 346
            self.match(bttestParser.T__38)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class TerminatorContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def NEWLINE(self):
            return self.getToken(bttestParser.NEWLINE, 0)

        def EOF(self):
            return self.getToken(bttestParser.EOF, 0)

        def getRuleIndex(self):
            return bttestParser.RULE_terminator

        def enterRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.enterTerminator(self)

        def exitRule(self, listener:ParseTreeListener):
            if isinstance( listener, bttestListener ):
                listener.exitTerminator(self)




    def terminator(self):

        localctx = bttestParser.TerminatorContext(self, self._ctx, self.state)
        self.enterRule(localctx, 46, self.RULE_terminator)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 348
            _la = self._input.LA(1)
            if not(_la==bttestParser.EOF or _la==bttestParser.NEWLINE):
                self._errHandler.recoverInline(self)
            else:
                self.consume()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx



    def sempred(self, localctx:RuleContext, ruleIndex:int, predIndex:int):
        if self._predicates == None:
            self._predicates = dict()
        self._predicates[17] = self.expression_sempred
        pred = self._predicates.get(ruleIndex, None)
        if pred is None:
            raise Exception("No predicate with index:" + str(ruleIndex))
        else:
            return pred(localctx, predIndex)

    def expression_sempred(self, localctx:ExpressionContext, predIndex:int):
            if predIndex == 0:
                return self.precpred(self._ctx, 7)
         

            if predIndex == 1:
                return self.precpred(self._ctx, 6)
         

            if predIndex == 2:
                return self.precpred(self._ctx, 5)
         

            if predIndex == 3:
                return self.precpred(self._ctx, 4)
         

            if predIndex == 4:
                return self.precpred(self._ctx, 3)
         

            if predIndex == 5:
                return self.precpred(self._ctx, 2)
         

            if predIndex == 6:
                return self.precpred(self._ctx, 1)
         

            if predIndex == 7:
                return self.precpred(self._ctx, 11)
         

            if predIndex == 8:
                return self.precpred(self._ctx, 10)
         




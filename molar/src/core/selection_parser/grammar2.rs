//##############################
//#  Selection grammar - alternative
//##############################

use crate::core::Vector3f;

struct AstNode {
    id: NodeId,
    payload: Payload,
}

impl AstNode {
    fn new(id: NodeId, payload: Payload) -> Self {
        Self {id, payload}
    }
}

enum Payload {
    Children(Vec<AstNode>),
    Str(String),
    Regex(regex::bytes::Regex),
    Int(i32),
    IntPair(i32,i32),
    Uint(usize),
    Float(f32),
    Char(char),
    Vec3(Vector3f),
    None,
}

enum NodeId {
    FnName,
    FnResname,
    FnResid,
    FnResindex,
    FnIndex,
    FnChain,

    FnModified,

    ArgStr,
    ArgRegex,
    ArgChar,
    ArgInt,
    ArgIntRange,
    ArgIntLt,
    ArgIntGt,
    ArgIntLeq,
    ArgIntGeq,

    Or,
    And,
    Not,

    SameResidue,
    SameChain,
    SameMolecule,
    Around,

    CmpEq,
    CmpNeq,
    CmpLt,
    CmpLeq,
    CmpGt,
    CmpGeq,
    CmpLtLt,
    CmpLeqLt,
    CmpLtLeq,
    CmpLeqLeq,
    CmpGtGt,
    CmpGeqGt,
    CmpGtGeq,
    CmpGeqGeq,

    Plus,
    Minus,
    UnaryMinus,
    Mul,
    Div,
    Pow,
    Rem,
    Abs,

    Vec3,
}

peg::parser! {
    pub(super) grammar selection_parser() for str {
        // Optional whitespace
        rule _ = (" " / "\t")*
        // Mandatory whitespace
        rule __ = (" " / "\t")+
        // Mandatory whitespace unless followed by paren
        rule ___ = _ &"(" / __

        rule uint() -> u32
            = n:$(['0'..='9']+)
            { n.parse().unwrap() }

        rule int()  -> i32
            = n:$(("-"/"+")? uint())
            { n.parse().unwrap() }

        rule float_val() -> f32
            = n:$((int() ("." uint())? / ("-"/"+") "." uint()) (("e"/"E") int())?)
            { n.parse().unwrap() }


        //================================================================
        
        // functions - rules that evaluate to list of indexes directly
        rule function() -> AstNode 
        = function_str() / function_int() / function_char()

        // Function names
        rule fn_name()      -> NodeId = "name"      {NodeId::FnName}
        rule fn_resname()   -> NodeId = "resname"   {NodeId::FnResname}
        rule fn_resid()     -> NodeId = "resid"     {NodeId::FnResid}
        rule fn_resindex()  -> NodeId = "resindex"  {NodeId::FnResindex}
        rule fn_index()     -> NodeId = "index"     {NodeId::FnIndex}
        rule fn_chain()     -> NodeId = "chain"     {NodeId::FnChain}

        // String functions
        rule function_str() -> AstNode
        = id:(fn_name() / fn_resname()) "(" args:function_args_str() ")"
        {
            AstNode::new(id, Payload::Children(args))
        }

        rule function_args_str() -> Vec<AstNode> 
        = (str_arg() / regex_arg()) ++ (__ / _ "," _);
        
        rule regex_arg() -> AstNode
        = "/" s:$((!"/" [_])+) "/"
        {?
            match regex::bytes::Regex::new(&format!("^{s}$")) {
                Ok(r) => Ok(
                    AstNode::new(NodeId::ArgRegex,Payload::Regex(r))
                ),
                Err(_) => Err("Invalid regex value"),
            }
        }

        rule str_arg() -> AstNode
        = s:$((![' '|'/'|')'] [_])+)
        { 
            AstNode::new(NodeId::ArgStr,Payload::Str(s.to_owned()))
        }

        // Int functions
        rule function_int() -> AstNode
        = id:(fn_resid() / fn_resindex() / fn_index()) "(" args:function_args_int() ")"
        {
            AstNode::new(id, Payload::Children(args))
        }

        rule function_args_int() -> Vec<AstNode> 
        = (int_range_arg() / int_arg() 
        / int_lt_arg() / int_leq_arg() 
        / int_gt_arg() / int_geq_arg()) ++ (__ / _ "," _);

        rule int_range_arg() -> AstNode
            = i1:int() _ ":" _ i2:int()
            { AstNode::new(NodeId::ArgIntRange, Payload::IntPair(i1,i2)) }

        rule int_arg() -> AstNode
            = i:int()
            { AstNode::new(NodeId::ArgInt, Payload::Int(i)) }

        rule int_lt_arg() -> AstNode
            = "<" _ i:int()
            { AstNode::new(NodeId::ArgIntLt, Payload::Int(i)) }

        rule int_leq_arg() -> AstNode
            = "<=" _ i:int()
            { AstNode::new(NodeId::ArgIntLeq, Payload::Int(i)) }

        rule int_gt_arg() -> AstNode
            = ">" _ i:int()
            { AstNode::new(NodeId::ArgIntGt, Payload::Int(i)) }

        rule int_geq_arg() -> AstNode
            = ">=" _ i:int()
            { AstNode::new(NodeId::ArgIntGeq, Payload::Int(i)) }

        // Character functions
        rule function_char() -> AstNode
        = id:fn_chain() "(" args:function_args_char() ")"
        {
            AstNode::new(id, Payload::Children(args))
        }

        rule function_args_char() -> Vec<AstNode> 
        = char_arg() ++ (__ / _ "," _);

        rule char_arg() -> AstNode
            = c:['a'..='z' | 'A'..='Z' | '0'..='9']
            { AstNode::new(NodeId::ArgChar, Payload::Char(c)) }
    
        // Modifiers - "methods" of functions called with ".method()", could be chained

        rule modifier() -> AstNode
        = same_modifier() / around_modifier();

        rule same_modifier() -> AstNode
        = "same(" _ id:(same_residue() / same_chain() / same_molecule()) _ ")" 
        { AstNode::new(id, Payload::None) }

        rule same_residue() -> NodeId
        = "residue" { NodeId::SameResidue }

        rule same_chain() -> NodeId
        = "chain" { NodeId::SameChain }

        rule same_molecule() -> NodeId
        = "molecule" { NodeId::SameMolecule }

        rule around_modifier() -> AstNode
        = "around(" d:float_val() ")"
        { AstNode::new(NodeId::Around, Payload::Float(d)) }

        // Comparisons


        // Logic
        pub rule logical() -> AstNode
        = precedence!{
            // Binary
            x:(@) (_ "||" _ / _ "or" ___) y:@ { AstNode::new(NodeId::Or, Payload::Children(vec![x,y])) }
            x:(@) (_ "&&" _ / _ "and" ___) y:@ { AstNode::new(NodeId::And, Payload::Children(vec![x,y])) }
            
            // Unary prefixes
            ("!" _ / "not" ___) v:@ { AstNode::new(NodeId::Not, Payload::Children(vec![v])) }
            x:@ "." m:modifier() { AstNode::new(NodeId::FnModified, Payload::Children(vec![x,m])) }
            //t:same_expr() ___ v:@ { LogicalNode::Same(t,Box::new(v)) }
            // Within from inner selection
            //p:within_expr() ___ v:@ {LogicalNode::Within(p,Box::new(v))}
            --
            //v:modifier_chain() {v}
            v:function() {v}
            //v:comparison_expr() { LogicalNode::Comparison(v) }
            //v:comparison_expr_chained() { LogicalNode::Comparison(v) }
            //v:compound() {LogicalNode::Compound(v)}
            // Within from point
            //p:within_expr() ___ v:vec3()  {LogicalNode::WithinPoint(p,v)}
            //"all" _ { LogicalNode::All }
            "(" _ e:logical() _ ")" { e }
        }
    } // grammar
} // parser


#[cfg(test)]
mod tests {
    use super::selection_parser::{logical};

    #[test]
    fn test1() {
        let _ast = logical("(resid(1 2 5:6 >=100) || !name(CA).same(residue).around(2.5)).same(molecule)")
            .expect("Error generating AST");
    }
}
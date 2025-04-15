//##############################
//#  Selection grammar - alternative
//##############################

use rayon::vec;

use crate::core::{PbcDims, Vector3f, PBC_FULL, PBC_NONE};

struct AstNode {
    id: NodeId,
    children: Vec<AstNode>,
}

impl AstNode {
    fn new(id: NodeId, children: Vec<AstNode>) -> Self {
        Self {id, children}
    }
}

enum NodeId {
    FnName,
    FnResname,
    FnResid,
    FnResindex,
    FnIndex,
    FnChain,

    FnModified,

    ArgStr(String),
    ArgRegex(regex::bytes::Regex),
    ArgChar(char),
    ArgInt(i32),
    ArgIntRange(i32,i32),
    ArgIntLt(i32),
    ArgIntGt(i32),
    ArgIntLeq(i32),
    ArgIntGeq(i32),

    Or,
    And,
    Not,

    SameResidue,
    SameChain,
    SameMolecule,
    Around(f32,PbcDims,bool),

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

    Float(f32),
    Vec3(Vector3f),
    Vdw,
    Mass,
    Charge,
    X,
    Y,
    Z,
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

        rule float() -> f32
            = n:$((int() ("." uint())? / ("-"/"+") "." uint()) (("e"/"E") int())?)
            { n.parse().unwrap() }


        //================================================================
        
        // functions - rules that evaluate to list of indexes
        rule function() -> AstNode 
        = function_str() / function_int() / function_char()

        // Function names
        rule fn_name()      -> NodeId = "name"      {NodeId::FnName}
        rule fn_resname()   -> NodeId = "resname"   {NodeId::FnResname}
        rule fn_resid()     -> NodeId = "resid"     {NodeId::FnResid}
        rule fn_resindex()  -> NodeId = "resindex"  {NodeId::FnResindex}
        rule fn_index()     -> NodeId = "index"     {NodeId::FnIndex}
        rule fn_chain()     -> NodeId = "chain"     {NodeId::FnChain}

        rule arg_delim()
        = __ / _ "," _ ;

        // String functions
        rule function_str() -> AstNode
        = id:(fn_name() / fn_resname()) "(" args:function_args_str() ")"
        {
            AstNode::new(id, args)
        }

        rule function_args_str() -> Vec<AstNode> 
        = (str_arg() / regex_arg()) ++ arg_delim();
        
        rule regex_arg() -> AstNode
        = "/" s:$((!"/" [_])+) "/"
        {?
            match regex::bytes::Regex::new(&format!("^{s}$")) {
                Ok(r) => Ok(
                    AstNode::new(NodeId::ArgRegex(r),vec![])
                ),
                Err(_) => Err("Invalid regex value"),
            }
        }

        rule str_arg() -> AstNode
        = s:$((![' '|'/'|')'] [_])+)
        { 
            AstNode::new(NodeId::ArgStr(s.to_owned()),vec![])
        }

        // Int functions
        rule function_int() -> AstNode
        = id:(fn_resid() / fn_resindex() / fn_index()) "(" args:function_args_int() ")"
        {
            AstNode::new(id, args)
        }

        rule function_args_int() -> Vec<AstNode> 
        = (int_range_arg() / int_arg() 
        / int_lt_arg() / int_leq_arg() 
        / int_gt_arg() / int_geq_arg()) ++ arg_delim();

        rule int_range_arg() -> AstNode
            = i1:int() _ ":" _ i2:int()
            { AstNode::new(NodeId::ArgIntRange(i1,i2),vec![]) }

        rule int_arg() -> AstNode
            = i:int()
            { AstNode::new(NodeId::ArgInt(i), vec![]) }

        rule int_lt_arg() -> AstNode
            = "<" _ i:int()
            { AstNode::new(NodeId::ArgIntLt(i), vec![]) }

        rule int_leq_arg() -> AstNode
            = "<=" _ i:int()
            { AstNode::new(NodeId::ArgIntLeq(i), vec![]) }

        rule int_gt_arg() -> AstNode
            = ">" _ i:int()
            { AstNode::new(NodeId::ArgIntGt(i), vec![]) }

        rule int_geq_arg() -> AstNode
            = ">=" _ i:int()
            { AstNode::new(NodeId::ArgIntGeq(i), vec![]) }

        // Character functions
        rule function_char() -> AstNode
        = id:fn_chain() "(" args:function_args_char() ")"
        {
            AstNode::new(id, args)
        }

        rule function_args_char() -> Vec<AstNode> 
        = char_arg() ++ arg_delim();

        rule char_arg() -> AstNode
            = c:['a'..='z' | 'A'..='Z' | '0'..='9']
            { AstNode::new(NodeId::ArgChar(c), vec![]) }
    
        // Modifiers - "methods" of functions called with ".method()", could be chained

        rule modifier() -> AstNode
        = same_modifier() / around_modifier();

        rule same_modifier() -> AstNode
        = "same(" _ id:(same_residue() / same_chain() / same_molecule()) _ ")" 
        { AstNode::new(id, vec![]) }

        rule same_residue() -> NodeId
        = "residue" { NodeId::SameResidue }

        rule same_chain() -> NodeId
        = "chain" { NodeId::SameChain }

        rule same_molecule() -> NodeId
        = "molecule" { NodeId::SameMolecule }

        rule around_modifier() -> AstNode
        = "around(" d:float() pbc:arg_pbc()? slf:arg_self()? ")"
        {
            let pbc = pbc.unwrap_or(PBC_NONE);
            let slf = slf.unwrap_or(true);
            AstNode::new(NodeId::Around(d,pbc,slf), vec![])
        }

        rule arg_pbc() -> PbcDims
        = arg_delim() p:pbc_dim()*<3> { PbcDims::new(p[0],p[1],p[2]) }

        rule pbc_dim() -> bool
        = v:$(['1'|'y'|'0'|'n']) {
            match v {
                "1" | "y" => true,
                "0" | "n" => false,
                _ => unreachable!()
            }
        }

        rule arg_self() -> bool
        = arg_delim() s:(self_literal() / noself_literal()) 
        { s }

        rule self_literal() -> bool
        = "self" { true }

        rule noself_literal() -> bool
        = "noself" { false }

        // Comparisons - expressions evaluating to the list indexes
        rule cmp_op_eq() ->  NodeId = "==" {NodeId::CmpEq}
        rule cmp_op_neq() -> NodeId = "!=" {NodeId::CmpNeq}
        rule cmp_op_leq() -> NodeId = "<=" {NodeId::CmpLeq}
        rule cmp_op_lt() ->  NodeId = "<"  {NodeId::CmpLt}
        rule cmp_op_geq() -> NodeId = ">=" {NodeId::CmpGeq}
        rule cmp_op_gt() ->  NodeId = ">"  {NodeId::CmpGt}

        // Simple comparison
        rule cmp_expr() -> AstNode =
            a:math_expr() _
            op:(cmp_op_eq() /cmp_op_neq()/
                cmp_op_leq()/cmp_op_lt() /
                cmp_op_geq()/cmp_op_gt() ) _
            b:math_expr() 
            {
                AstNode::new(op, vec![a,b])    
            }
        
        // Math expressions on float or Vec3 arguments
        rule math_expr() -> AstNode
        = precedence!{
            x:(@) _ "+" _ y:@ { AstNode::new(NodeId::Plus,vec![x,y]) }
            x:(@) _ "-" _ y:@ { AstNode::new(NodeId::Minus,vec![x,y]) }
            --
            x:(@) _ "*" _ y:@ { AstNode::new(NodeId::Mul,vec![x,y]) }
            x:(@) _ "/" _ y:@ { AstNode::new(NodeId::Div,vec![x,y]) }
            --
            x:@ _ "^" _ y:(@) { AstNode::new(NodeId::Pow,vec![x,y]) }
            --
            "-" _ v:@ { AstNode::new(NodeId::UnaryMinus,vec![v]) }
            --
            v:float() { AstNode::new(NodeId::Float(v), vec![]) }
            ['x'|'X'] { AstNode::new(NodeId::X,        vec![]) }
            ['y'|'Y'] { AstNode::new(NodeId::Y,        vec![]) }
            ['z'|'Z'] { AstNode::new(NodeId::Z,        vec![]) }
            "vdw"     { AstNode::new(NodeId::Vdw,      vec![]) }
            "mass"    { AstNode::new(NodeId::Mass,     vec![]) }
            "charge"  { AstNode::new(NodeId::Charge,   vec![]) }
            // d:distance() {MathNode::Dist(d)}
            // keyword_occupancy() { MathNode::Occupancy }
            // keyword_bfactor() { MathNode::Bfactor }
            // f:math_function_name() _ "(" _ e:math_expr() _ ")" {
            //     MathNode::Function(f,Box::new(e))
            // }
            "(" _ e:math_expr() _ ")" { e }
        }
        

        // Logic
        pub rule logical() -> AstNode
        = precedence!{
            // Binary
            x:(@) (_ "||" _ / _ "or" ___) y:@ { AstNode::new(NodeId::Or, vec![x,y]) }
            x:(@) (_ "&&" _ / _ "and" ___) y:@ { AstNode::new(NodeId::And, vec![x,y]) }
            
            // Unary prefixes
            ("!" _ / "not" ___) v:@ { AstNode::new(NodeId::Not, vec![v]) }
            x:@ "." m:modifier() { AstNode::new(NodeId::FnModified, vec![x,m]) }
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
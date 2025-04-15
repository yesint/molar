//##############################
//#  Selection grammar - alternative
//##############################

use crate::core::{PbcDims, Vector3f, PBC_NONE};

#[derive(Clone)]
struct AstNode {
    id: NodeId,
    children: Vec<AstNode>,
}

impl AstNode {
    fn new(id: impl Into<NodeId>, children: Vec<AstNode>) -> Self {
        Self { id: id.into(), children }
    }
}

#[derive(Clone)]
enum NodeId {
    Index(IndexNodeId),
    Float(FloatNodeId),
    Vec(VecNodeId),
}

#[derive(Clone)]
enum IndexNodeId {
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
    ArgIntRange(i32, i32),
    ArgIntLt(i32),
    ArgIntGt(i32),
    ArgIntLeq(i32),
    ArgIntGeq(i32),

    Or,
    And,
    Not,

    Around(f32, PbcDims, bool),
    AroundPoint,

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

    SameResidue,
    SameChain,
    SameMolecule,
}

#[derive(Clone)]
enum FloatNodeId {
    FoldFloat,
    
    Plus,
    Minus,
    UnaryMinus,
    Mul,
    Div,
    Pow,
    Rem,
    Abs,

    Float(f32),
    Vdw,
    Mass,
    Charge,
    X,
    Y,
    Z,
    VecElement(usize),

    Line,
    Plane,
}

#[derive(Clone)]
enum VecNodeId {
    FoldVec3,
    Com(PbcDims),
    Cog(PbcDims),
    Vec3(Vector3f),
}

impl From<IndexNodeId> for NodeId {
    fn from(value: IndexNodeId) -> Self {
        Self::Index(value)
    }
}

impl From<FloatNodeId> for NodeId {
    fn from(value: FloatNodeId) -> Self {
        Self::Float(value)
    }
}

impl From<VecNodeId> for NodeId {
    fn from(value: VecNodeId) -> Self {
        Self::Vec(value)
    }
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
        rule fn_name()      -> IndexNodeId = "name"      {IndexNodeId::FnName}
        rule fn_resname()   -> IndexNodeId = "resname"   {IndexNodeId::FnResname}
        rule fn_resid()     -> IndexNodeId = "resid"     {IndexNodeId::FnResid}
        rule fn_resindex()  -> IndexNodeId = "resindex"  {IndexNodeId::FnResindex}
        rule fn_index()     -> IndexNodeId = "index"     {IndexNodeId::FnIndex}
        rule fn_chain()     -> IndexNodeId = "chain"     {IndexNodeId::FnChain}

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
                    AstNode::new(IndexNodeId::ArgRegex(r),vec![])
                ),
                Err(_) => Err("Invalid regex value"),
            }
        }

        rule str_arg() -> AstNode
        = s:$((![' '|'/'|')'] [_])+)
        {
            AstNode::new(IndexNodeId::ArgStr(s.to_owned()),vec![])
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
            { AstNode::new(IndexNodeId::ArgIntRange(i1,i2),vec![]) }

        rule int_arg() -> AstNode
            = i:int()
            { AstNode::new(IndexNodeId::ArgInt(i), vec![]) }

        rule int_lt_arg() -> AstNode
            = "<" _ i:int()
            { AstNode::new(IndexNodeId::ArgIntLt(i), vec![]) }

        rule int_leq_arg() -> AstNode
            = "<=" _ i:int()
            { AstNode::new(IndexNodeId::ArgIntLeq(i), vec![]) }

        rule int_gt_arg() -> AstNode
            = ">" _ i:int()
            { AstNode::new(IndexNodeId::ArgIntGt(i), vec![]) }

        rule int_geq_arg() -> AstNode
            = ">=" _ i:int()
            { AstNode::new(IndexNodeId::ArgIntGeq(i), vec![]) }

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
            { AstNode::new(IndexNodeId::ArgChar(c), vec![]) }

        // Modifiers - "methods" of functions called with ".method()", could be chained

        rule modifier() -> AstNode
        = same_modifier() / around_modifier();

        rule same_modifier() -> AstNode
        = "same(" _ id:(same_residue() / same_chain() / same_molecule()) _ ")"
        { AstNode::new(id, vec![]) }

        rule same_residue() -> IndexNodeId
        = "residue" { IndexNodeId::SameResidue }

        rule same_chain() -> IndexNodeId
        = "chain" { IndexNodeId::SameChain }

        rule same_molecule() -> IndexNodeId
        = "molecule" { IndexNodeId::SameMolecule }

        rule around_modifier() -> AstNode
        = "around(" d:float() pbc:arg_pbc()? slf:arg_self()? ")"
        {
            let pbc = pbc.unwrap_or(PBC_NONE);
            let slf = slf.unwrap_or(true);
            AstNode::new(IndexNodeId::Around(d,pbc,slf), vec![])
        }

        // Aggregators - compute a float or vector from indexes
        rule fold_to_vec3() -> AstNode
        = com() / cog()

        rule com() -> AstNode
        = "com(" pbc:pbc_dims()? ")"
        {
            let pbc = pbc.unwrap_or(PBC_NONE);
            AstNode::new(VecNodeId::Com(pbc), vec![])
        }

        rule cog() -> AstNode
        = "cog(" pbc:pbc_dims()? ")"
        {
            let pbc = pbc.unwrap_or(PBC_NONE);
            AstNode::new(VecNodeId::Cog(pbc), vec![])
        }

        rule arg_pbc() -> PbcDims
        = arg_delim() p:pbc_dims() { p }

        rule pbc_dims() -> PbcDims
        = p:pbc_dim()*<3> { PbcDims::new(p[0],p[1],p[2]) }

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

        // Vectors
        #[cache_left_rec]
        rule vec3() -> AstNode
        = vec3_lit() / vec3_folded() 

        rule vec3_lit() -> AstNode
        = "[" v:(float()**<3> arg_delim()) "]"
        {
            AstNode::new(VecNodeId::Vec3(Vector3f::new(v[0],v[1],v[2])), vec![])
        }

        rule line() -> AstNode
        = "line(" p1:vec3() arg_delim() p2:vec3() ")"
        { AstNode::new(FloatNodeId::Line, vec![p1,p2]) }

        #[cache_left_rec]
        rule vec3_folded() -> AstNode
        = l:logical() "." f:fold_to_vec3()
        { AstNode::new(VecNodeId::FoldVec3, vec![l,f]) }

        // Comparisons - expressions evaluating to the list indexes
        rule cmp_op_eq() ->  IndexNodeId = "==" {IndexNodeId::CmpEq}
        rule cmp_op_neq() -> IndexNodeId = "!=" {IndexNodeId::CmpNeq}
        rule cmp_op_leq() -> IndexNodeId = "<=" {IndexNodeId::CmpLeq}
        rule cmp_op_lt() ->  IndexNodeId = "<"  {IndexNodeId::CmpLt}
        rule cmp_op_geq() -> IndexNodeId = ">=" {IndexNodeId::CmpGeq}
        rule cmp_op_gt() ->  IndexNodeId = ">"  {IndexNodeId::CmpGt}

        // Simple comparison
        #[cache_left_rec]
        rule cmp_expr() -> AstNode =
            a:math_expr() _
            op:(cmp_op_eq() /cmp_op_neq()/
                cmp_op_leq()/cmp_op_lt() /
                cmp_op_geq()/cmp_op_gt() ) _
            b:math_expr()
            {
                AstNode::new(op, vec![a,b])
            }

        // Math expressions on float arguments
        #[cache_left_rec]
        rule math_expr() -> AstNode
        = precedence!{
            x:(@) _ "+" _ y:@ { AstNode::new(FloatNodeId::Plus,vec![x,y]) }
            x:(@) _ "-" _ y:@ { AstNode::new(FloatNodeId::Minus,vec![x,y]) }
            --
            x:(@) _ "*" _ y:@ { AstNode::new(FloatNodeId::Mul,vec![x,y]) }
            x:(@) _ "/" _ y:@ { AstNode::new(FloatNodeId::Div,vec![x,y]) }
            --
            x:@ _ "^" _ y:(@) { AstNode::new(FloatNodeId::Pow,vec![x,y]) }
            --
            "-" _ v:@ { AstNode::new(FloatNodeId::UnaryMinus,vec![v]) }
            --
            v:float() { AstNode::new(FloatNodeId::Float(v), vec![]) }
            ['x'|'X'] { AstNode::new(FloatNodeId::X,        vec![]) }
            ['y'|'Y'] { AstNode::new(FloatNodeId::Y,        vec![]) }
            ['z'|'Z'] { AstNode::new(FloatNodeId::Z,        vec![]) }
            "vdw"     { AstNode::new(FloatNodeId::Vdw,      vec![]) }
            "mass"    { AstNode::new(FloatNodeId::Mass,     vec![]) }
            "charge"  { AstNode::new(FloatNodeId::Charge,   vec![]) }
            v:vec3() ".x" {  AstNode::new(FloatNodeId::VecElement(0),vec![v])  }
            v:vec3() ".y" {  AstNode::new(FloatNodeId::VecElement(1),vec![v])  }
            v:vec3() ".z" {  AstNode::new(FloatNodeId::VecElement(2),vec![v])  }

            "(" _ e:math_expr() _ ")" { e }
        }

        // Logic
        #[cache_left_rec]
        pub rule logical() -> AstNode
        = precedence!{
            // Binary
            x:(@) (_ "||" _ / _ "or" ___) y:@ { AstNode::new(IndexNodeId::Or, vec![x,y]) }
            x:(@) (_ "&&" _ / _ "and" ___) y:@ { AstNode::new(IndexNodeId::And, vec![x,y]) }

            // Unary prefixes
            ("!" _ / "not" ___) v:@ { AstNode::new(IndexNodeId::Not, vec![v]) }

            x:@ "." f:fold_to_vec3() "." a:around_modifier() {
                AstNode::new(
                    IndexNodeId::AroundPoint,
                    vec![AstNode::new(VecNodeId::FoldVec3, vec![x,f]), a]
                )
            }

            x:@ "." m:modifier() { AstNode::new(IndexNodeId::FnModified, vec![x,m]) }
            --

            v:function() { v }
            v:cmp_expr() { v }

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
    use super::selection_parser::logical;

    // resid[1,2 5:6,>=100] || !name[CA CB].same(residue).around(2.5,yyy,self)).same(molecule)

    #[test]
    fn test1() {
        let _ast =
            logical("(resid(1 2 5:6 >=100) || !name(CA).same(residue).around(2.5)).same(molecule)")
                .expect("Error generating AST");
    }

    #[test]
    fn test2() {
        let _ast = logical("name(A).com().around(2.5).same(chain)").expect("Error generating AST");
    }

    #[test]
    fn test3() {
        let _ast = logical("((x < name(A).com().x + 2) || name(A)).around(2.5)").expect("Error generating AST");
    }
}

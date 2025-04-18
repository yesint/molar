//##############################
//#  Selection grammar - alternative
//##############################

use std::{cell::RefCell, collections::HashMap};

use crate::core::{PbcDims, Vector3f, PBC_NONE};

#[derive(Clone,Debug)]
struct AstNode {
    id: NodeId,
    children: Vec<AstNode>,
}

fn new_node(id: impl Into<NodeId>, children: Vec<AstNode>) -> AstNode {
    AstNode { id: id.into(), children }
}

#[derive(Clone,Debug)]
enum NodeId {
    // Nodes evaluating to array of indexes
    Index(IndexNodeId),
    // Nodes evaluating to float
    Float(FloatNodeId),
    // Nodes evaluating to Vec3
    Vec(VecNodeId),
    // Nodes that don't evaluate by they own
    Arg(ArgNodeId),
    // Variables
    Var(VarNodeId),
}

#[derive(Clone,Debug)]
enum IndexNodeId {
    TopLevel,

    FnName,
    FnResname,
    FnResid,
    FnResindex,
    FnIndex,
    FnChain,
    FnModified,

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

#[derive(Clone,Debug)]
enum ArgNodeId {
    ArgStr(String),
    ArgRegex(regex::bytes::Regex),
    ArgChar(char),
    ArgInt(i32),
    ArgIntRange(i32, i32),
    ArgIntLt(i32),
    ArgIntGt(i32),
    ArgIntLeq(i32),
    ArgIntGeq(i32),
}

#[derive(Clone,Debug)]
enum VarNodeId {
    Logical(String),
    Float(String),
    Vec(String),
}

#[derive(Clone,Debug)]
enum FloatNodeId {
    FoldFloat,
    
    Add,
    Sub,
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

    Func,
    FuncAbs,
    FuncSqrt,
}

#[derive(Clone,Debug)]
enum VecNodeId {
    FoldVec3,
    Com(PbcDims),
    Cog(PbcDims),
    Vec3(Vector3f),

    Add,
    Sub,
    UnaryMinus,
    MulScalar,
    DivScalar,
    PowScalar,
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

impl From<ArgNodeId> for NodeId {
    fn from(value: ArgNodeId) -> Self {
        Self::Arg(value)
    }
}

impl From<VarNodeId> for NodeId {
    fn from(value: VarNodeId) -> Self {
        Self::Var(value)
    }
}

peg::parser! {
    pub(super) grammar selection_parser(symbol_table: &RefCell<HashMap<String,AstNode>>) for str {
        use IndexNodeId as I;
        use FloatNodeId as F;
        use VecNodeId as V;
        use ArgNodeId as A;

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
        rule fn_name()      -> I = "name"      {I::FnName}
        rule fn_resname()   -> I = "resname"   {I::FnResname}
        rule fn_resid()     -> I = "resid"     {I::FnResid}
        rule fn_resindex()  -> I = "resindex"  {I::FnResindex}
        rule fn_index()     -> I = "index"     {I::FnIndex}
        rule fn_chain()     -> I = "chain"     {I::FnChain}

        rule arg_delim()
        = __ / _ "," _ ;

        // String functions
        rule function_str() -> AstNode
        = id:(fn_name() / fn_resname()) "(" args:function_args_str() ")"
        {
            new_node(id, args)
        }

        rule function_args_str() -> Vec<AstNode>
        = (str_arg() / regex_arg()) ++ arg_delim();

        rule regex_arg() -> AstNode
        = "/" s:$((!"/" [_])+) "/"
        {?
            match regex::bytes::Regex::new(&format!("^{s}$")) {
                Ok(r) => Ok(
                    new_node(A::ArgRegex(r),vec![])
                ),
                Err(_) => Err("Invalid regex"),
            }
        }

        rule str_arg() -> AstNode
        = s:$((![' '|'/'|')'] [_])+)
        {
            new_node(A::ArgStr(s.to_owned()),vec![])
        }

        // Int functions
        rule function_int() -> AstNode
        = id:(fn_resid() / fn_resindex() / fn_index()) "(" args:function_args_int() ")"
        {
            new_node(id, args)
        }

        rule function_args_int() -> Vec<AstNode>
        = (int_range_arg() / int_arg()
        / int_lt_arg() / int_leq_arg()
        / int_gt_arg() / int_geq_arg()) ++ arg_delim();

        rule int_range_arg() -> AstNode
            = i1:int() _ ":" _ i2:int()
            { new_node(A::ArgIntRange(i1,i2),vec![]) }

        rule int_arg() -> AstNode
            = i:int()
            { new_node(A::ArgInt(i), vec![]) }

        rule int_lt_arg() -> AstNode
            = "<" _ i:int()
            { new_node(A::ArgIntLt(i), vec![]) }

        rule int_leq_arg() -> AstNode
            = "<=" _ i:int()
            { new_node(A::ArgIntLeq(i), vec![]) }

        rule int_gt_arg() -> AstNode
            = ">" _ i:int()
            { new_node(A::ArgIntGt(i), vec![]) }

        rule int_geq_arg() -> AstNode
            = ">=" _ i:int()
            { new_node(A::ArgIntGeq(i), vec![]) }

        // Character functions
        rule function_char() -> AstNode
        = id:fn_chain() "(" args:function_args_char() ")"
        {
            new_node(id, args)
        }

        rule function_args_char() -> Vec<AstNode>
        = char_arg() ++ arg_delim();

        rule char_arg() -> AstNode
            = c:['a'..='z' | 'A'..='Z' | '0'..='9']
            { new_node(A::ArgChar(c), vec![]) }

        // Modifiers - "methods" of functions called with ".method()", could be chained

        rule modifier() -> AstNode
        = same_modifier() / around_modifier();

        rule same_modifier() -> AstNode
        = "same(" _ id:(same_residue() / same_chain() / same_molecule()) _ ")"
        { new_node(id, vec![]) }

        rule same_residue() -> I
        = "residue" { I::SameResidue }

        rule same_chain() -> I
        = "chain" { I::SameChain }

        rule same_molecule() -> I
        = "molecule" { I::SameMolecule }

        rule around_modifier() -> AstNode
        = "around(" d:float() pbc:arg_pbc()? slf:arg_self()? ")"
        {
            let pbc = pbc.unwrap_or(PBC_NONE);
            let slf = slf.unwrap_or(true);
            new_node(I::Around(d,pbc,slf), vec![])
        }

        // Compute a vector from index rule
        rule fold_to_vec3() -> AstNode
        = com() / cog()

        rule com() -> AstNode
        = "com(" pbc:pbc_dims()? ")"
        {
            let pbc = pbc.unwrap_or(PBC_NONE);
            new_node(V::Com(pbc), vec![])
        }

        rule cog() -> AstNode
        = "cog(" pbc:pbc_dims()? ")"
        {
            let pbc = pbc.unwrap_or(PBC_NONE);
            new_node(V::Cog(pbc), vec![])
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

        // Literal vector: [1.0 2.0 3.0]
        rule vec3_lit() -> AstNode
        = "[" v:(float()**<3> arg_delim()) "]"
        {
            new_node(V::Vec3(Vector3f::new(v[0],v[1],v[2])), vec![])
        }

        rule line() -> AstNode
        = "line(" p1:vec3() arg_delim() p2:vec3() ")"
        { new_node(F::Line, vec![p1,p2]) }

        // Vector obtained by "folding" an index rule: (...).com()
        #[cache_left_rec]
        rule vec3_folded() -> AstNode
        = l:logical() "." f:fold_to_vec3()
        { new_node(V::FoldVec3, vec![l,f]) }

        // Comparisons - expressions evaluating to the list indexes
        rule cmp_op_eq()  -> I = "==" {I::CmpEq}
        rule cmp_op_neq() -> I = "!=" {I::CmpNeq}
        rule cmp_op_leq() -> I = "<=" {I::CmpLeq}
        rule cmp_op_lt()  -> I = "<"  {I::CmpLt}
        rule cmp_op_geq() -> I = ">=" {I::CmpGeq}
        rule cmp_op_gt()  -> I = ">"  {I::CmpGt}

        // Simple comparison
        #[cache_left_rec]
        rule cmp_expr() -> AstNode =
            a:math_expr() _
            op:(cmp_op_eq() /cmp_op_neq()/
                cmp_op_leq()/cmp_op_lt() /
                cmp_op_geq()/cmp_op_gt() ) _
            b:math_expr()
            {
                new_node(op, vec![a,b])
            }


        rule math_function() -> AstNode
        = math_func_abs() / math_func_sqrt()

        rule math_func_abs() -> AstNode
        = "abs(" _ ")"
        { new_node(F::FuncAbs, vec![]) }

        rule math_func_sqrt() -> AstNode
        = "abs(" _ ")"
        { new_node(F::FuncAbs, vec![]) }

        // Math expressions on float arguments
        #[cache_left_rec]
        rule math_expr() -> AstNode
        = precedence!{
            x:(@) _ "+" _ y:@ { new_node(F::Add, vec![x,y]) }
            x:(@) _ "-" _ y:@ { new_node(F::Sub, vec![x,y]) }
            --
            x:(@) _ "*" _ y:@ { new_node(F::Mul, vec![x,y]) }
            x:(@) _ "/" _ y:@ { new_node(F::Div, vec![x,y]) }
            --
            x:@ _ "^" _ y:(@) { new_node(F::Pow, vec![x,y]) }
            --
            // Math functions
            x:@ "." f:math_function() { new_node(F::Func, vec![x,f]) }
            --
            "-" _ v:@ { new_node(F::UnaryMinus,vec![v]) }
            --
            v:float()     { new_node(F::Float(v), vec![]) }
            ['x'|'X']     { new_node(F::X,        vec![]) }
            ['y'|'Y']     { new_node(F::Y,        vec![]) }
            ['z'|'Z']     { new_node(F::Z,        vec![]) }
            "vdw"         { new_node(F::Vdw,      vec![]) }
            "mass"        { new_node(F::Mass,     vec![]) }
            "charge"      { new_node(F::Charge,   vec![]) }
            v:vector_math_expr() ".x" { new_node(F::VecElement(0),vec![v]) }
            v:vector_math_expr() ".y" { new_node(F::VecElement(1),vec![v]) }
            v:vector_math_expr() ".z" { new_node(F::VecElement(2),vec![v]) }
            v:var_float() { v }

            "(" _ e:math_expr() _ ")" { e }
        }

        #[cache_left_rec]
        rule vector_math_expr() -> AstNode
        = precedence!{
            x:(@) _ "+" _ y:@ { new_node(V::Add, vec![x,y]) }
            x:(@) _ "-" _ y:@ { new_node(V::Sub, vec![x,y]) }
            --
            x:@ _ "*" _ s:math_expr() { new_node(V::MulScalar, vec![x,s]) }
            x:@ _ "/" _ s:math_expr() { new_node(V::DivScalar, vec![x,s]) }
            --
            x:@ _ "^" _ s:math_expr() { new_node(V::PowScalar, vec![x,s]) }
            --
            "-" _ v:@ { new_node(V::UnaryMinus,vec![v]) }
            --
            v:vec3() { v }
            v:var_vec() { v }
            l:logical() "." f:fold_to_vec3() {
                new_node(V::FoldVec3, vec![l,f])
            }

            "(" _ e:vector_math_expr() _ ")" { e }
        }

        // Variables
        rule var_index() -> AstNode
        = "$" i:ident()
        {?
            if let Some(var) = symbol_table.borrow().get(&i) {
                if let NodeId::Index(_) = var.id {
                    Ok(var.clone())
                } else {
                    Err("variable is not a logical expression")    
                }
            } else {
                Err("undeclared variable")
            }
        }

        rule var_float() -> AstNode
        = "$" i:ident()
        {?
            if let Some(var) = symbol_table.borrow().get(&i) {
                if let NodeId::Float(_) = var.id {
                    Ok(var.clone())
                } else {
                    Err("variable is not a float expression")    
                }
            } else {
                Err("undeclared variable")
            }
        }

        rule var_vec() -> AstNode
        = "$" i:ident()
        {?
            if let Some(var) = symbol_table.borrow().get(&i) {
                if let NodeId::Vec(_) = var.id {
                    Ok(var.clone())
                } else {
                    Err("variable is not a vector expression")    
                }
            } else {
                Err("undeclared variable")
            }
        }

        rule var_def()
        = "$" i:ident() _ "=" _ e:(vector_math_expr() / math_expr() / logical()) _ ";"
        { symbol_table.borrow_mut().insert(i,e); }

        rule ident() -> String 
        = s:$(['a'..='z'|'A'..='Z'|'_'] ['a'..='z'|'A'..='Z'|'_'|'0'..='9']*) {s.to_owned()}

        // Logic
        #[cache_left_rec]
        rule logical() -> AstNode
        = precedence!{
            // Binary
            x:(@) (_ "||" _ / _ "or" ___)  y:@ { new_node(I::Or, vec![x,y]) }
            x:(@) (_ "&&" _ / _ "and" ___) y:@ { new_node(I::And, vec![x,y]) }
            ("!" _ / "not" ___) v:@ { new_node(I::Not, vec![v]) }
            --
            x:@ "." m:modifier() { new_node(I::FnModified, vec![x,m]) }
            //--
            v:function() { v }
            v:cmp_expr() { v }
            v:var_index() { v }
            v:vector_math_expr() "." a:around_modifier() {
                new_node(I::AroundPoint, vec![v,a])
            }
            "(" _ e:logical() _ ")" { e }
        }

        pub rule top_level() -> AstNode
        = (var_def() ** _) _ l:logical()
        {
            new_node(I::TopLevel, vec![l])
        }

    } // grammar
} // parser

#[cfg(test)]
mod tests {
    use std::{cell::RefCell, collections::HashMap};

    use super::selection_parser::top_level;

    // resid[1,2 5:6,>=100] || !name[CA CB].same(residue).around(2.5,yyy,self)).same(molecule)

    #[test]
    fn test1() {
        let table = RefCell::new(HashMap::new());
        let _ast =
        // !name(CA /N.*/).same(residue).around(2.5 yyy noself).same(molecule)
        // not same molecule as within 2.5 of same residue as name CA
            top_level("(resid(1 2 5:6 >=100) || !name(CA).same(residue).around(2.5)).same(molecule)", &table)
                .expect("Error generating AST");
    }

    #[test]
    fn test2() {
        let table = RefCell::new(HashMap::new());
        let _ast = top_level("name(A).com().around(2.5).same(chain)", &table).expect("Error generating AST");
    }

    #[test]
    fn test3() {
        let table = RefCell::new(HashMap::new());
        let _ast = top_level("((x < name(A).com().x + 2) || name(A)).around(2.5)",&table).unwrap();
    }

    #[test]
    fn test4() {
        let table = RefCell::new(HashMap::new());
        let _ast = top_level("([1 2 3]*7+name(CA).com()^2)*5.around(2.5)", &table).unwrap();
    }

    #[test]
    fn test_vars() {
        let table = RefCell::new(HashMap::new());
        let ast = top_level("$ca=resname(ALA); name(CB) && ([1 1 1]+$ca.com()).x>0", &table).unwrap();
        println!("{:?}",ast);
    }

    #[test]
    fn pass_vars() {
        let table = RefCell::new(HashMap::new());
        let ast = top_level("$ca=resname(ALA); name(CB) && ([1 1 1]+$ca.com()).x>0", &table).unwrap();
        println!("{:?}",ast);
    }
}

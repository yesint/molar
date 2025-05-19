//##############################
//#  Selection grammar - alternative
//##############################

use std::{cell::RefCell, collections::HashMap};

use crate::core::{PbcDims, Vector3f, PBC_NONE};


peg::parser! {
    pub(super) grammar selection_parser() for str {
        //================================================
        // Terminals
        //================================================

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
        
        rule arg_delim()
            = __ / _ "," _ ;    

        //================================================
        // Top level
        //================================================

        pub rule top_level()
        = index_expr()
        
        //================================================
        // Index
        //================================================
    
        #[cache_left_rec]
        rule index_expr()
        = precedence!{
            // Binary
            x:(@) (_ "||" _ / _ "or" ___)  y:@ { }
            x:(@) (_ "&&" _ / _ "and" ___) y:@ { }
            --
            ("!" _ / "not" ___) v:@ { }
            --
            vec_expr() "." around() { }
            index_expr() "." around() { }
            // Special case of immediate .com().around()
            index_expr() "." center() "." around() { }
            
            function() { }
            "(" _  index_expr() _ ")" { }
        }

        //================================================
        // Vector
        //================================================

        #[cache_left_rec]
        rule vec_expr()
        = precedence!{
            x:(@) _ "+" _ y:@ { }
            x:(@) _ "-" _ y:@ { }
            --
            "[" _ float() _ float() _ float() _ "]" { }
            
            "(" _  vec_expr() _ ")" { }
            index_expr() ".com()" { }
        }

        //================================================
        // Float expr
        //================================================

        rule float_expr()
        = float()

        //================================================
        // Functions
        //================================================

        rule function()
        = function_str() / function_int() / function_char()

        // Function names
        // rule fn_name()      -> I = "name"      {I::FnName}
        // rule fn_resname()   -> I = "resname"   {I::FnResname}
        // rule fn_resid()     -> I = "resid"     {I::FnResid}
        // rule fn_resindex()  -> I = "resindex"  {I::FnResindex}
        // rule fn_index()     -> I = "index"     {I::FnIndex}
        // rule fn_chain()     -> I = "chain"     {I::FnChain}

        rule fn_name()      = "name"    
        rule fn_resname()   = "resname" 
        rule fn_resid()     = "resid"   
        rule fn_resindex()  = "resindex"
        rule fn_index()     = "index"   
        rule fn_chain()     = "chain"   

        // String functions
        rule function_str()
        = (fn_name() / fn_resname()) "(" function_args_str() ")"
        // {
        //     new_node(id, args)
        // }

        rule function_args_str()
        = (str_arg() / regex_arg()) ++ arg_delim();

        rule regex_arg()
        = "/" s:$((!"/" [_])+) "/"
        // {?
        //     match regex::bytes::Regex::new(&format!("^{s}$")) {
        //         Ok(r) => Ok(
        //             new_node(A::ArgRegex(r),vec![])
        //         ),
        //         Err(_) => Err("Invalid regex"),
        //     }
        // }

        rule str_arg()
        = ((![' '|'/'|')'] [_])+)
        // {
        //     new_node(A::ArgStr(s.to_owned()),vec![])
        // }

        // Int functions
        rule function_int()
        = (fn_resid() / fn_resindex() / fn_index()) "(" function_args_int() ")"
        // {
        //     new_node(id, args)
        // }

        rule function_args_int()
        = (int_range_arg() / int_arg()
        / int_lt_arg() / int_leq_arg()
        / int_gt_arg() / int_geq_arg()) ++ arg_delim();

        rule int_range_arg()
            = i1:int() _ ":" _ i2:int()
            //{ new_node(A::ArgIntRange(i1,i2),vec![]) }

        rule int_arg()
            = i:int()
            //{ new_node(A::ArgInt(i), vec![]) }

        rule int_lt_arg()
            = "<" _ i:int()
            //{ new_node(A::ArgIntLt(i), vec![]) }

        rule int_leq_arg()
            = "<=" _ i:int()
            //{ new_node(A::ArgIntLeq(i), vec![]) }

        rule int_gt_arg()
            = ">" _ i:int()
            //{ new_node(A::ArgIntGt(i), vec![]) }

        rule int_geq_arg()
            = ">=" _ i:int()
            //{ new_node(A::ArgIntGeq(i), vec![]) }

        // Character functions
        rule function_char()
        = fn_chain() "(" function_args_char() ")"
        // {
        //     new_node(id, args)
        // }

        rule function_args_char()
        = char_arg() ++ arg_delim();

        rule char_arg()
            = c:['a'..='z' | 'A'..='Z' | '0'..='9']
            //{ new_node(A::ArgChar(c), vec![]) }

        //================================================
        // Around
        //================================================

        rule around()
        = "around(" d:float() pbc:arg_pbc()? slf:arg_self()? ")" { }

        rule center()
        = com() / cog()

        rule com()
        = "com(" pbc_dims()? ")"
        // {
        //     let pbc = pbc.unwrap_or(PBC_NONE);
        //     new_node(V::Com(pbc), vec![])
        // }

        rule cog()
        = "cog(" pbc_dims()? ")"
        // {
        //     let pbc = pbc.unwrap_or(PBC_NONE);
        //     new_node(V::Cog(pbc), vec![])
        // }

        //================================================
        // PBC
        //================================================

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


    } // grammar
} // parser

#[cfg(test)]
mod tests {
    use std::{cell::RefCell, collections::HashMap};

    use super::selection_parser::top_level;

    // resid[1,2 5:6,>=100] || !name[CA CB].same(residue).around(2.5,yyy,self)).same(molecule)

    #[test]
    fn test1() {
        let _ast =
        // !name(CA /N.*/).same(residue).around(2.5 yyy noself).same(molecule)
        // not same molecule as within 2.5 of same residue as name CA
            top_level("((name(CA)||name(CB)).com()+name(CC).com()).around(3)").unwrap();
            //top_level("name(CA).com().around()").unwrap();
    }

    
}

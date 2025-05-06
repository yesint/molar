//##############################
//#  Selection grammar
//##############################

use super::ast::*;
use crate::prelude::*;

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

        rule float() -> MathNode
            = v:float_val() { MathNode::Float(v) }

        // Keywords
        rule keyword_name() -> AtomAttr = "name" {AtomAttr::Name}
        rule keyword_resid() -> AtomAttr = "resid" {AtomAttr::Resid}
        rule keyword_resindex() -> AtomAttr = "resindex" {AtomAttr::Resindex}
        rule keyword_resname() -> AtomAttr = "resname" {AtomAttr::Resname}
        rule keyword_index() -> AtomAttr = "index" {AtomAttr::Index}
        rule keyword_chain() -> AtomAttr = "chain" {AtomAttr::Chain}
        rule keyword_residue() -> AtomAttr = "residue" {AtomAttr::Residue}
        rule keyword_occupancy() -> AtomAttr = ("occupancy" / "occ") {AtomAttr::Occupancy}
        rule keyword_bfactor() -> AtomAttr = ("bfactor" / "beta") {AtomAttr::Bfactor}

        rule int_keyword_expr() -> KeywordNode
            = s:(keyword_resid() / keyword_resindex() / keyword_index()) __
              v:((int_range() / int_single()) ++ __)
            {
                match s {
                    AtomAttr::Resid => KeywordNode::Resid(v),
                    AtomAttr::Resindex => KeywordNode::Resindex(v),
                    AtomAttr::Index => KeywordNode::Index(v),
                    _ => unreachable!()
                }
            }

        rule int_range() -> IntKeywordArg
            = i1:int() _ ":" _ i2:int()
            { IntKeywordArg::IntRange(i1,i2) }

        rule int_single() -> IntKeywordArg
            = i:int()
            { IntKeywordArg::Int(i) }

        rule chain_keyword_expr() -> KeywordNode
        = s:keyword_chain() __ v:(['a'..='z' | 'A'..='Z' | '0'..='9'] ++ __)
        {
            KeywordNode::Chain(v)
        }

        rule str_keyword_expr() -> KeywordNode
        = s:(keyword_name() / keyword_resname()) __
            v:((str_value() / regex_value()) ++ __)
        {?
            match s {
                AtomAttr::Name => Ok(KeywordNode::Name(v)),
                AtomAttr::Resname => Ok(KeywordNode::Resname(v)),
                _ => Err(unreachable!())
            }
        }

        rule keyword_expr() -> KeywordNode
        = v:(int_keyword_expr() / str_keyword_expr() / chain_keyword_expr()) {v}

        rule regex_value() -> StrKeywordArg
        = "/" s:$((!"/" [_])+) "/"
        {?
            match regex::bytes::Regex::new(&format!("^{s}$")) {
                Ok(r) => Ok(StrKeywordArg::Regex(r)),
                Err(_) => Err("Invalid regex value"),
            }
        }

        rule str_value() -> StrKeywordArg
        = !("and"/"or") s:$((![' '|'/'|')'] [_])+)
        { StrKeywordArg::Str(s.to_owned()) }


        // 3-vector value
        rule vec3() -> VectorNode = vec3_spaces() / vec3_comas() / vec3_com() / vec3_cog()

        rule vec3_spaces() -> VectorNode
        = x:float_val() __ y:float_val() __ z:float_val() {
            VectorNode::Const(Pos::new(x, y, z))
        }

        rule vec3_comas() -> VectorNode
        = "[" _ x:float_val() _ "," _ y:float_val() _ "," _ z:float_val() _ "]" {
            VectorNode::Const(Pos::new(x, y, z))
        }

        rule vec3_com() -> VectorNode
        = "com" __ p:pbc_expr()? "of" ___ v:logical_expr() {
            let pbc = match p {
                Some(dims) => dims,
                None => PBC_NONE,
            };
            VectorNode::Com(v.into(),pbc)
        }

        rule vec3_cog() -> VectorNode
        = "com" __ p:pbc_expr()? "of" ___ v:logical_expr() {
            let pbc = match p {
                Some(dims) => dims,
                None => PBC_NONE,
            };
            VectorNode::Cog(v.into(),pbc)
        }

        //rule nth_of() -> Pos
        //= logical_expr()

        // Distance
        rule distance() -> DistanceNode
        = distance_point() / distance_line_2points() / distance_line_point_dir()
        / distance_plane_3points() / distance_plane_point_normal()

        rule distance_point() -> DistanceNode
        = "dist" __ b:pbc_expr()? "point" __ p:vec3() {
            DistanceNode::Point(p,b.unwrap_or(PBC_NONE))
        }

        rule distance_line_2points() -> DistanceNode
        = "dist" __ b:pbc_expr()? "line" __ p1:vec3() __ p2:vec3() {
            DistanceNode::Line(p1.into(),p2.into(),b.unwrap_or(PBC_NONE))
        }

        rule distance_line_point_dir() -> DistanceNode
        = "dist" __ b:pbc_expr()? "line" __ p:vec3() __ "dir" __ dir:vec3() {
            DistanceNode::LineDir(p.into(),dir.into(),b.unwrap_or(PBC_NONE))
        }

        rule distance_plane_3points() -> DistanceNode
        = "dist" __ b:pbc_expr()? "plane" __ p1:vec3() __ p2:vec3() __ p3:vec3() {
            DistanceNode::Plane(p1.into(),p2.into(),p3.into(),b.unwrap_or(PBC_NONE))
        }

        rule distance_plane_point_normal() -> DistanceNode
        = "dist" __ b:pbc_expr()? "plane" __ p:vec3() __ "normal" __ n:vec3() {
            DistanceNode::PlaneNormal(p.into(),n.into(),b.unwrap_or(PBC_NONE))
        }

        rule abs_function() -> MathFunctionName = "abs" {MathFunctionName::Abs}
        rule sqrt_function() -> MathFunctionName = "sqrt" {MathFunctionName::Sqrt}
        rule sin_function() -> MathFunctionName = "sin" {MathFunctionName::Sin}
        rule cos_function() -> MathFunctionName = "cos" {MathFunctionName::Cos}

        rule math_function_name() -> MathFunctionName
        = abs_function() / sqrt_function() / sin_function() / cos_function()

        // Math
        rule math_expr() -> MathNode
        = precedence!{
            x:(@) _ "+" _ y:@ { MathNode::Add(x.into(), y.into()) }
            x:(@) _ "-" _ y:@ { MathNode::Sub(x.into(), y.into()) }
                    "-" _ v:@ { MathNode::Neg(v.into()) }
            --
            x:(@) _ "*" _ y:@ { MathNode::Mul(x.into(), y.into()) }
            x:(@) _ "/" _ y:@ { MathNode::Div(x.into(), y.into()) }
            --
            x:@ _ "^" _ y:(@) { MathNode::Pow(x.into(), y.into()) }
            --
            v:float() {v}
            ['x'|'X'] { MathNode::X }
            ['y'|'Y'] { MathNode::Y }
            ['z'|'Z'] { MathNode::Z }
            "vdw" {MathNode::Vdw}
            "mass" {MathNode::Mass}
            "charge" {MathNode::Charge}
            d:distance() {MathNode::Dist(d)}
            keyword_occupancy() { MathNode::Occupancy }
            keyword_bfactor() { MathNode::Bfactor }
            f:math_function_name() _ "(" _ e:math_expr() _ ")" {
                MathNode::Function(f, e.into())
            }
            "(" _ e:math_expr() _ ")" { e }
        }

        // Comparisons
        rule comparison_op_eq() -> ComparisonOp = "==" {ComparisonOp::Eq}
        rule comparison_op_neq() -> ComparisonOp = "!=" {ComparisonOp::Neq}
        rule comparison_op_leq() -> ComparisonOp = "<=" {ComparisonOp::Leq}
        rule comparison_op_lt() -> ComparisonOp = "<" {ComparisonOp::Lt}
        rule comparison_op_geq() -> ComparisonOp = ">=" {ComparisonOp::Geq}
        rule comparison_op_gt() -> ComparisonOp = ">" {ComparisonOp::Gt}

        // Simple comparison
        rule comparison_expr() -> ComparisonNode =
            a:math_expr() _
            op:(comparison_op_eq()/comparison_op_neq()/
                comparison_op_leq()/comparison_op_lt()/
                comparison_op_geq()/comparison_op_gt()) _
            b:math_expr() {
                use ComparisonOp as C;
                match op {
                    C::Eq => { ComparisonNode::Eq(a,b) },
                    C::Neq => { ComparisonNode::Neq(a,b) },
                    C::Leq => { ComparisonNode::Leq(a,b) },
                    C::Lt => { ComparisonNode::Lt(a,b) },
                    C::Geq => { ComparisonNode::Geq(a,b) },
                    C::Gt => { ComparisonNode::Gt(a,b) },
                    _ => unreachable!(),
                }
            }

        // Chained comparison
        rule comparison_expr_chained() -> ComparisonNode
        = comparison_expr_chained_l() / comparison_expr_chained_r()

        rule comparison_expr_chained_l() -> ComparisonNode
        =   a:math_expr() _
            op1:(comparison_op_leq()/comparison_op_lt()) _
            b:math_expr() _
            op2:(comparison_op_leq()/comparison_op_lt()) _
            c:math_expr()
        {
            use ComparisonOp as C;
            match (op1,op2) {
                (C::Lt,C::Lt) => { ComparisonNode::LtLt(a,b,c) },
                (C::Lt,C::Leq) => { ComparisonNode::LtLeq(a,b,c) },
                (C::Leq,C::Lt) => { ComparisonNode::LeqLt(a,b,c) },
                (C::Leq,C::Leq) => { ComparisonNode::LeqLeq(a,b,c) },
                _ => unreachable!(),
            }
        }

        rule comparison_expr_chained_r() -> ComparisonNode
        =   a:math_expr() _
            op1:(comparison_op_geq()/comparison_op_gt()) _
            b:math_expr() _
            op2:(comparison_op_geq()/comparison_op_gt()) _
            c:math_expr()
        {
            use ComparisonOp as C;
            match (op1,op2) {
                (C::Gt,C::Gt) => { ComparisonNode::GtGt(a,b,c) },
                (C::Gt,C::Geq) => { ComparisonNode::GtGeq(a,b,c) },
                (C::Geq,C::Gt) => { ComparisonNode::GeqGt(a,b,c) },
                (C::Geq,C::Geq) => { ComparisonNode::GeqGeq(a,b,c) },
                _ => unreachable!(),
            }
        }

        // "Same" expressions
        rule same_expr() -> SameAttr
        = "same" __ t:(keyword_residue() / keyword_chain()) __ "as" {
            match t {
                AtomAttr::Residue => SameAttr::Residue,
                AtomAttr::Chain => SameAttr::Chain,
                _ => unreachable!(),
            }
        }

        // Single PBC dimention
        rule pbc_dim() -> bool
        = v:$("1" / "0" / "y" / "n") {
            match v {
                "1" | "y" => true,
                "0" | "n" => false,
                _ => unreachable!()
            }
        }

        // PBC for within
        rule pbc_expr() -> PbcDims
        = pbc_with_dims() / pbc_full_no_dims() / pbc_none_no_dims()

        rule pbc_full_no_dims() -> PbcDims
        = "pbc" __
        {PBC_FULL}

        rule pbc_none_no_dims() -> PbcDims
        = "nopbc" __
        {PBC_NONE}

        rule pbc_with_dims() -> PbcDims
        = "pbc" __ p:(pbc_dim()*<3>) __ {
            PbcDims::new(p[0],p[1],p[2])
        }

        // Within
        rule within_expr() -> WithinParams
        = "within" __ d:float_val() __ p:pbc_expr()? s:$(("self" __)?) "of" {
            let pbc = match p {
                Some(dims) => dims,
                None => PBC_NONE,
            };
            let include_inner = !s.is_empty();
            WithinParams {cutoff: d, pbc, include_inner}
        }

        // COM
        rule com_expr() -> PbcDims
        = "com" __ p:pbc_expr()? "of" {
            match p {
                Some(dims) => dims,
                None => PBC_NONE,
            }
        }

        // COG
        rule cog_expr() -> PbcDims
        = ("cog" / "center") __ p:pbc_expr()? "of" {
            match p {
                Some(dims) => dims,
                None => PBC_NONE,
            }
        }

        pub rule compound() -> ChemicalNode 
        = protein() / backbone() / sidechain() / 
          water() / not_water() / hydrogen() / not_hydrogen ();

        pub rule protein() -> ChemicalNode = "protein" _ { ChemicalNode::Protein };
        pub rule sidechain() -> ChemicalNode = "sidechain" _ { ChemicalNode::Sidechain };
        pub rule backbone() -> ChemicalNode = "backbone" _ { ChemicalNode::Backbone };
        pub rule water() -> ChemicalNode = "water" _ { ChemicalNode::Water };
        pub rule not_water() -> ChemicalNode = "now" _ { ChemicalNode::NotWater };
        pub rule hydrogen() -> ChemicalNode = "hydrogen" _ { ChemicalNode::Hydrogen };
        pub rule not_hydrogen() -> ChemicalNode = "noh" _ { ChemicalNode::NotHydrogen };

        // Logic
        pub rule logical_expr() -> LogicalNode
        = precedence!{
            // Binary
            x:(@) _ "or" _ y:@ { LogicalNode::Or(x.into(), y.into()) }
            x:(@) _ "and" _ y:@ { LogicalNode::And(x.into(), y.into()) }
            // Unary prefixes
            "not" ___ v:@ { LogicalNode::Not(v.into()) }
            t:same_expr() ___ v:@ { LogicalNode::Same(t, v.into()) }
            // Within from inner selection
            p:within_expr() ___ v:@ {LogicalNode::Within(p, v.into())}
            --
            v:keyword_expr() { LogicalNode::Keyword(v) }
            v:comparison_expr() { LogicalNode::Comparison(v) }
            v:comparison_expr_chained() { LogicalNode::Comparison(v) }
            v:compound() {LogicalNode::Chemical(v)}
            // Within from point
            p:within_expr() ___ v:vec3()  {LogicalNode::WithinPoint(p,v)}
            "all" _ { LogicalNode::All }
            "(" _ e:logical_expr() _ ")" { e }
        }
    } // grammar
} // parser

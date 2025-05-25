use std::{collections::HashMap, path::Path};
use itertools::Itertools;

use crate::prelude::*;

/// Representation of Gromacs index files
pub struct NdxFile {
    groups: HashMap<String, SVec>,
}

impl NdxFile {
    /// Creates a new NdxFile by parsing a Gromacs index file
    pub fn new(path: impl AsRef<Path>) -> Result<Self, SelectionError> {
        let path = path.as_ref();
        let ndx_str = std::fs::read_to_string(path)
            .map_err(|e| NdxError::NdxIo(path.to_owned(), e))?;
        
        let mut groups = HashMap::new();
        let mut current_group = None;
        let mut current_numbers = Vec::new();

        for line in ndx_str.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            if line.starts_with('[') && line.ends_with(']') {
                // Store previous group if exists
                if let Some(group_name) = current_group.take() {
                    if !current_numbers.is_empty() {
                        groups.insert(group_name, current_numbers.into());
                    }
                }
                // Start new group
                let group_name = line[1..line.len()-1].trim().to_string();
                current_group = Some(group_name);
                current_numbers = Vec::new();
            } else if let Some(group_name) = &current_group {
                // Parse numbers for current group
                current_numbers.extend(
                    line.split_whitespace()
                        .map(|s| s.parse::<usize>())
                        .map_ok(|i| i - 1) // Convert to zero-based
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(|e| NdxError::Parse(group_name.clone(), e))?,
                );
            } else {
                return Err(NdxError::MalformedNdxFile(path.to_owned()))?;
            }
        }

        // Store last group if exists
        if let Some(group_name) = current_group {
            if !current_numbers.is_empty() {
                groups.insert(group_name, current_numbers.into());
            }
        }

        Ok(Self { groups })
    }

    /// Get an index group by name
    pub fn get_group(&self, name: impl AsRef<str>) -> Result<&SVec, SelectionError> {
        let gr = name.as_ref();
        Ok(self.groups.get(gr).ok_or_else(|| NdxError::NoGroup(gr.to_owned()))?)
    }
}

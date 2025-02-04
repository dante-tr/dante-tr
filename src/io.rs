use crate::TandemRepeat;
use crate::Module;

pub fn get_modules(
    left_flank: &[u8], repeat: &TandemRepeat, right_flank: &[u8]
) -> Vec<Module> {
    let mut modules = Vec::new();
    modules.push(left_flank.into());
    modules_add_motif(&mut modules, repeat);
    modules.push(right_flank.into());
    return modules;
}

fn modules_add_motif(modules: &mut Vec<Module>, motif: &TandemRepeat) {
    for i in 0..motif.copy_unit.len() {
        modules.push((&motif.copy_unit[i][..], motif.copy_number[i]).into())
    }
}


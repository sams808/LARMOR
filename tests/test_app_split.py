"""The G5 contract: ``MainWindow`` is a facade over the ``mw_*`` mixins.

``larmor/desktop/app.py`` keeps ``__init__``, the two Qt event overrides and
``main()``; every other method lives in one ``larmor/desktop/mw_<block>.py``
mixin and the QThreads in ``larmor/desktop/workers.py``. These tests pin what
the split must preserve and what pyflakes cannot see:

- the frozen member surface (every method and class constant stays reachable
  on the class, so ``win.<name>`` and ``MainWindow.<staticmethod>`` keep
  working and a quiet fall-through to a Qt attribute of the same name is
  caught),
- the facade holding nothing but construction and the two Qt overrides,
- mixins that are plain Python -- no ``__init__``, no ``Signal``, no Qt event
  override -- with pairwise-disjoint names (a method left in two modules
  would be shadowed silently by MRO order),
- every ``self.<name>(`` call resolving to a method or to an attribute that
  is assigned on ``self`` somewhere (the "pyflakes for self" the gate lacks;
  it covers menu slots no test clicks),
- no mixin importing ``larmor.desktop.app`` (app imports the mixins at module
  level, so the reverse would be circular),
- the names other modules import from ``larmor.desktop.app``,
- the menu tree, the shortcut set and the toolbar action count of an
  offscreen window.

When a menu label or a shortcut changes ON PURPOSE, regenerate GOLDEN_MENU
and SHORTCUTS: the failure message prints the fresh rows ready to paste.
"""
import ast
import os
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[1]
DESKTOP = ROOT / "larmor" / "desktop"

#: what app.py itself may define on MainWindow
FACADE_METHODS = frozenset({"__init__", "keyPressEvent", "closeEvent"})
#: Qt virtuals a mixin must never define (base order decides which wins)
QT_OVERRIDES = frozenset({
    "keyPressEvent", "keyReleaseEvent", "closeEvent", "resizeEvent",
    "showEvent", "hideEvent", "event", "eventFilter", "mousePressEvent",
    "mouseReleaseEvent", "mouseMoveEvent", "paintEvent", "dragEnterEvent",
    "dropEvent", "changeEvent", "focusInEvent", "focusOutEvent"})
#: names other code imports from larmor.desktop.app (cli, launcher, three
#: dialogs' lazy ``_load_any``, the tests) and where each is defined now
REEXPORTS = ("MainWindow", "main", "asset_path", "TUTORIALS", "FitWorker",
             "Fit2DWorker", "KernelWarmWorker", "SimWorker", "_emit_progress",
             "humanize_error", "_fit_tol", "_load_any")
HOMES = {"FitWorker": "larmor.desktop.workers",
         "Fit2DWorker": "larmor.desktop.workers",
         "KernelWarmWorker": "larmor.desktop.workers",
         "SimWorker": "larmor.desktop.workers",
         "_emit_progress": "larmor.desktop.workers",
         "humanize_error": "larmor.desktop.workers",
         "_fit_tol": "larmor.desktop.workers",
         "_load_any": "larmor.desktop.mw_files"}
#: children of these submenus come from QSettings (recent files / recipes),
#: which LARMOR_NO_SESSION does not gate, so they vary per machine
EXCLUDE_CHILDREN = frozenset({"Open &recent", "&Apply recipe"})

# 317 MainWindow members
SURFACE = (
    'remove_sites', 'add_sidebands_for_line',
    '__init__', '_build_menus', '_menu', '_add', '_MODEL_GROUPS',
    '_rebuild_recent', '_add_recent',
    '_add_recent_recipe', '_rebuild_apply_recipe', 'apply_recipe_browse',
    '_recipe_model_from', 'show_chi2_map', 'compare_with_saved_fit',
    'apply_recipe', 'add_sidebands', '_toggle_watch', '_watch_targets',
    '_retarget_watch', '_watched_changed', '_reload_watched',
    '_active_plot_widget', 'copy_plot', 'save_plot_image', 'open_integrals',
    'open_nmr_table', 'open_convert', 'edit_processing_steps',
    'edit_computing_params', 'edit_mqmas_f1_ref', '_open_manual',
    'open_command_palette', '_open_tutorial', '_about', '_show_more',
    '_build_toolbar', '_build_sidebar', '_toggle_scroll_nudge',
    '_toggle_ref_ranges', '_update_ref_ranges', '_build_czjzek_display_menu',
    '_set_czjzek_display', '_build_explorer_dock', '_open_source_paths',
    '_on_explorer_renamed', '_build_workspaces_dock',
    '_build_datasets_dock', '_build_bottom_docks', '_build_right_dock',
    '_fit_to_screen', '_build_panels_menu', '_toggle_health_strip',
    '_update_enabled', '_build_axis_unit_menu', '_set_axis_unit',
    '_apply_axis_unit', '_format_x', '_update_exp_label', '_update_mas_label',
    'edit_fit_tol', 'save_constraint_set', 'apply_constraint_set',
    'restrict_glass_protocol', 'edit_experiment', '_baseline_mode',
    'start_twopoint_bg', '_twopoint_mode', 'apply_twopoint_bg',
    'apply_manual_baseline', 'apply_iterbaseline', 'add_zone', 'clear_zones',
    '_sync_zones', '_zones_changed', 'zoom_full', 'zoom_sites', 'autoscale_y',
    '_toggle_resid', '_toggle_comp', '_toggle_labels', '_toggle_paddles',
    '_build_theme_menu', '_set_aesthetic_override', '_build_textsize_menu',
    '_set_text_size', '_apply_theme_live', '_set_theme',
    '_maybe_show_welcome', 'open_file', 'open_expno', 'open_varian',
    'open_sample', 'open_fid', '_fid_to_workbench', '_fid_to_2d', '_last_dir',
    '_suggest_dir', '_doc_title', 'save_project', 'open_project',
    '_open_project_path', '_confirm_open_mode', '_report_project_notes',
    '_restore_1d_entry', '_restore_2d_entry', '_snapshot_doc', '_apply_doc',
    '_sync_active', '_register_ws', '_add_session_entry', '_refresh_ws_panel',
    'switch_workspace', 'close_workspace', 'open_workspace_entry',
    'save_workspace', '_resolve_open_path', '_offer_fxml_no_data',
    'load_source', '_load_source_body', '_load_nonfittable', '_update_sn',
    'add_overlay_dialog', '_read_overlay_source', '_add_overlay',
    'overlay_set_color', 'overlay_remove', 'overlay_visibility',
    'overlay_make_active', '_refresh_overlays', 'PROJ_COLOR',
    'overlay_1d_on_2d', '_current_1d_source', 'load_projection_1d',
    '_explorer_open', '_read_1d', 'back_to_2d', '_show_2d', '_warm_kernel',
    '_display_1d', '_trace_to_workbench', 'new_fit', '_recipe_copy',
    '_capture_state', '_restore_state', 'snapshot', 'undo', 'redo',
    '_set_add_mode', 'keyPressEvent', '_NUCLEUS_START',
    '_seed_nucleus_defaults', '_seed_from_data', '_remembered_site_defaults',
    '_remember_site_defaults', 'add_site_at', 'add_function_line',
    'autopick_lines', 'label_from_literature', 'predict_at_field',
    'add_background_spectrum', '_append_spectrum_site',
    'add_current_spectrum_line', '_SSB_SELF_MODELS', '_ssb_key',
    '_model_has_manifold', '_run_sideband_detection',
    '_maybe_offer_sidebands', 'detect_sidebands', '_show_sideband_offer',
    '_dismiss_sideband_offer', '_toggle_ssb_offer', '_fresh_params',
    '_ssb_set_rate_if_due', '_apply_sideband_choice',
    '_sideband_parent_index', '_add_sideband_manifold',
    '_add_sideband_copies', '_add_sideband_model_line', 'add_site_2d',
    'on_site_structure', '_move_site', '_remap_exprs_after_delete',
    'on_structure_changed', 'on_params_changed', '_update_paddles',
    'on_paddle_moved', 'on_paddle_released', '_kernel_progress_tick',
    'request_simulation', '_simulate_now', '_sim_busy_on', '_sim_busy_off',
    '_on_sim_failed', '_sim_done', 'run_fit', '_fit_frame', '_MODELS_2D',
    'run_fit_2d', '_fit2d_done', '_sanitize_constraints_before_fit',
    '_fit_failed', '_residual_noise_ratio', '_health_from_result',
    '_health_live', '_health_apply', '_health_show', '_health_reset',
    'show_fit_health', '_health_focus_param', '_health_show_residual',
    '_fit_done', 'run_quantify', 'copy_csv', 'copy_latex', 'copy_methods',
    'export_publication_bundle', 'apply_processing', 'start_phase_drag',
    '_phase_drag_mode', 'on_phase_dragged', 'on_phase_drag_released',
    '_on_pivot_moved', '_proc_spec_valid', '_mirror_view_actions',
    '_display_fallback', '_on_view_experiment_set', '_refresh_display',
    '_on_proc_view_changed', '_toggle_time_domain', '_set_channel',
    '_cycle_channel', 'start_calibrate', '_shift_recipe_positions',
    'on_calibrate_picked', 'toggle_measure', 'on_measure_changed',
    'reset_processing', 'save_recipe', '_guard_write', 'save_fit_as',
    'open_figure_dialog', 'open_wurst_correct', 'open_subtract',
    '_subtract_applied', 'save_spectrum', '_build_cofit_page',
    '_build_progress', '_show_fit_buttons', '_interrupt_fit',
    '_style_progress', '_progress_start', '_progress_tick', '_progress_end',
    '_COFIT_LABEL', '_cofit_tieable', '_cofit_tie_rebuild',
    '_cofit_tie_toggled', '_cofit_copy_param', '_cofit_home', 'open_cofit',
    '_default_tie', '_on_central_changed', '_cofit_split_even',
    '_set_cofit_dock', '_cofit_rebuild_tables', '_cofit_on_edit',
    '_cofit_struct', '_cofit_add_dataset', '_cofit_load',
    '_cofit_refresh_panels', '_cofit_active', '_cofit_simulate',
    '_cofit_simulate_now', '_cofit_sim_2d', 'run_cofit_fit',
    '_cofit_apply_to_main', 'close_cofit', 'open_per_site_relaxation',
    'open_vt', 'show_correlations', 'show_czjzek_dist',
    'open_qcpmg_batch_fields', 'open_qcpmg_fields', 'open_staticct',
    '_staticct_seed', 'open_herzfeld_berger', '_hb_seed',
    'open_referencing_audit', 'open_acquisition_table', 'open_session_inventory',
    '_apply_sr_correction', 'open_vocs',
    '_vocs_to_workbench', 'open_qcpmg', 'open_satrec', 'open_redor',
    'open_magres', '_magres_add_sites', 'open_twod', 'run_auto_fit',
    'run_errors_analysis',
    'run_monte_carlo', 'run_batch_report', 'run_batch_fit', 'run_seq_fit',
    'open_plotting_studio', 'plot_current_spectrum', '_open_session_entry',
    '_remember_figure', '_remember_batch', '_open_batch_session',
    '_session_file', '_persist_session', '_flush_session', '_autosave_tick',
    'closeEvent', '_restore_session', 'on_family_changed',
)

# 219 menu rows: (menu path, text, shortcut, checkable, has submenu)
GOLDEN_MENU = (
    ((), '&File', '', False, True),
    (('&File',), '&Open…', 'Ctrl+O', False, False),
    (('&File',), 'Open &sample…', 'Ctrl+Shift+S', False, False),
    (('&File',), 'Open &EXPNO / folder…', 'Ctrl+Shift+O', False, False),
    (('&File',), 'Open &FID…', 'Ctrl+F', False, False),
    (('&File',), 'Open &Varian / Agilent…', '', False, False),
    (('&File',), 'Open &recent', '', False, True),
    (('&File',), 'O&verlay a spectrum…', 'Ctrl+Shift+A', False, False),
    (('&File',), '&Watch the source file', '', True, False),
    (('&File',), 'Open pro&ject…', '', False, False),
    (('&File',), 'Save projec&t…', 'Ctrl+Alt+S', False, False),
    (('&File',), '&Save fit', 'Ctrl+S', False, False),
    (('&File',), 'Save fit &as…', 'Ctrl+Shift+E', False, False),
    (('&File',), 'E&xport', '', False, True),
    (('&File', 'E&xport'), 'Save s&pectrum as CSV…', '', False, False),
    (('&File', 'E&xport'), '&Figure…', '', False, False),
    (('&File', 'E&xport'), 'Save plot &image…', '', False, False),
    (('&File', 'E&xport'), '&Copy plot to clipboard', 'Ctrl+Shift+C', False, False),
    (('&File', 'E&xport'), 'Copy report &table as CSV', '', False, False),
    (('&File', 'E&xport'), 'Copy &LaTeX table', '', False, False),
    (('&File', 'E&xport'), '&Publication bundle…', '', False, False),
    (('&File',), '&Quit', '', False, False),
    ((), '&Edit', '', False, True),
    (('&Edit',), '↩  Undo', 'Ctrl+Z', False, False),
    (('&Edit',), '↪  Redo', 'Ctrl+Y', False, False),
    (('&Edit',), '&New fit', '', False, False),
    (('&Edit',), 'Add &line', '', False, True),
    (('&Edit', 'Add &line'), 'Gauss/Lorentz', '', True, False),
    (('&Edit', 'Add &line'), 'Gauss/Lorentz (area)', '', True, False),
    (('&Edit', 'Add &line'), 'Voigt (true)', '', True, False),
    (('&Edit', 'Add &line'), 'J-multiplet', '', True, False),
    (('&Edit', 'Add &line'), 'Quad CT (2nd order)', '', True, False),
    (('&Edit', 'Add &line'), 'Quad 1st order (satellites)', '', True, False),
    (('&Edit', 'Add &line'), 'Quad CT + CSA', '', True, False),
    (('&Edit', 'Add &line'), 'Czjzek (quad. distribution)', '', True, False),
    (('&Edit', 'Add &line'), 'ext. Czjzek', '', True, False),
    (('&Edit', 'Add &line'), 'Czjzek, general d  (GIM at d = 5)', '', True, False),
    (('&Edit', 'Add &line'), 'Czjzek + δiso–C_Q correlation', '', True, False),
    (('&Edit', 'Add &line'), 'Amorphous (Gaussian Cq/eta dist.)', '', True, False),
    (('&Edit', 'Add &line'), 'CSA powder (MAS/static)', '', True, False),
    (('&Edit', 'Add &line'), 'CSA distribution (disordered)', '', True, False),
    (('&Edit', 'Add &line'), 'Spinning sidebands', '', True, False),
    (('&Edit', 'Add &line'), 'Two-site exchange  (Bloch–McConnell)', '', True, False),
    (('&Edit', 'Add &line'), 'Function fit', '', True, False),
    (('&Edit',), 'Add a line at every &peak…', '', False, False),
    (('&Edit',), 'Add f&unction line…', '', False, False),
    (('&Edit',), 'Add background &spectrum…', '', False, False),
    (('&Edit',), '&Apply recipe', '', False, True),
    (('&Edit',), 'Spinning &sidebands', '', False, True),
    (('&Edit', 'Spinning &sidebands'), '&Detect spinning sidebands', 'Ctrl+Shift+D', False, False),
    (('&Edit', 'Spinning &sidebands'), 'Add spinning &sidebands…', '', False, False),
    (('&Edit', 'Spinning &sidebands'), 'Add a &copy of this spectrum…', '', False, False),
    (('&Edit', 'Spinning &sidebands'), '&Offer spinning-sideband detection on load', '', True, False),
    (('&Edit',), '&Constraints', '', False, True),
    (('&Edit', '&Constraints'), '&Label lines from literature ranges', '', False, False),
    (('&Edit', '&Constraints'), '&Restrict around current values…', '', False, False),
    (('&Edit', '&Constraints'), 'Save &constraints as…', '', False, False),
    (('&Edit', '&Constraints'), 'Apply saved co&nstraints…', '', False, False),
    (('&Edit',), 'Add fit &zone', '', False, False),
    (('&Edit',), 'Clear zones', '', False, False),
    ((), '&Process', '', False, True),
    (('&Process',), '&Experiment parameters…', '', False, False),
    (('&Process',), 'Processing s&teps…', '', False, False),
    (('&Process',), 'Show processing &panel', '', False, False),
    (('&Process',), '&Phase', '', False, True),
    (('&Process', '&Phase'), '&Autophase (ACME)', '', False, False),
    (('&Process', '&Phase'), '&Drag to phase', 'Ctrl+P', False, False),
    (('&Process',), '&Baseline', '', False, True),
    (('&Process', '&Baseline'), '&Polynomial (order 3)', '', False, False),
    (('&Process', '&Baseline'), '&Iterative (dead-time; Yon 2020)…', '', False, False),
    (('&Process', '&Baseline'), '&2-point background…', '', False, False),
    (('&Process', '&Baseline'), 'Subtract &averages', '', False, False),
    (('&Process',), '&Reference', '', False, True),
    (('&Process', '&Reference'), '&Calibrate axis…', '', False, False),
    (('&Process', '&Reference'), '&Measure Δ (ppm / Hz)', '', True, False),
    (('&Process', '&Reference'), 'Referencing a&udit…', '', False, False),
    (('&Process',), 'Re&gion / algebra', '', False, True),
    (('&Process', 'Re&gion / algebra'), '&Integrals && measurements…', '', False, False),
    (('&Process', 'Re&gion / algebra'), 'Subtract a spectrum (&background)…', '', False, False),
    (('&Process', 'Re&gion / algebra'), '&WURST excitation profile…', '', False, False),
    (('&Process', 'Re&gion / algebra'), 'Stitch frequency-stepped (&VOCS) spectra…', '', False, False),
    (('&Process',), '&FID ⇄ spectrum', 'Ctrl+T', True, False),
    (('&Process',), 'Display &channel', '', False, True),
    (('&Process', 'Display &channel'), '&Real', '', True, False),
    (('&Process', 'Display &channel'), '&Imaginary', '', True, False),
    (('&Process', 'Display &channel'), '&Magnitude |S|', '', True, False),
    (('&Process', 'Display &channel'), 'C&ycle channel', 'Ctrl+I', False, False),
    (('&Process',), 'Reset to &original', '', False, False),
    ((), 'F&it', '', False, True),
    (('F&it',), '&Simulate', 'F9', False, False),
    (('F&it',), '&Fit', 'F5', False, False),
    (('F&it',), '&Auto fit…', '', False, False),
    (('F&it',), '&Report', 'F6', False, False),
    (('F&it',), 'Fit &health details…', 'F7', False, False),
    (('F&it',), '&Errors', '', False, True),
    (('F&it', '&Errors'), 'χ² &profile…', '', False, False),
    (('F&it', '&Errors'), 'Monte-&Carlo errors…', '', False, False),
    (('F&it', '&Errors'), 'Parameter correlations…', '', False, False),
    (('F&it', '&Errors'), 'χ² map (parameter pair)…', '', False, False),
    (('F&it',), 'Compare with a saved fit…', '', False, False),
    (('F&it',), 'Co-&fit datasets…', '', False, False),
    (('F&it',), '&Predict at another field…', '', False, False),
    (('F&it',), '&MQMAS', '', False, True),
    (('F&it', '&MQMAS'), '2D MQMAS &viewer / fit…', '', False, False),
    (('F&it', '&MQMAS'), 'MQMAS F1 &reference…', '', False, False),
    (('F&it',), 'Fit se&ttings', '', False, True),
    (('F&it', 'Fit se&ttings'), 'Computing &parameters…', '', False, False),
    (('F&it', 'Fit se&ttings'), 'Fit completion &threshold…', '', False, False),
    (('F&it', 'Fit se&ttings'), '&Animate fits', '', True, False),
    ((), '&Series', '', False, True),
    (('&Series',), 'Batch &fit spectra…', '', False, False),
    (('&Series',), 'Se&quential fit…', '', False, False),
    (('&Series',), '&Batch fit report…', '', False, False),
    (('&Series',), '&Session inventory…', '', False, False),
    (('&Series',), 'E&xperimental section…', '', False, False),
    (('&Series',), '&Compare acquisition parameters…', '', False, False),
    ((), '&Tools', '', False, True),
    (('&Tools',), '&NMR table…', '', False, False),
    (('&Tools',), '&Conversion tools…', '', False, False),
    (('&Tools',), '&Herzfeld–Berger sideband analysis…', '', False, False),
    (('&Tools',), '&Read static pattern (C_Q, η)…', '', False, False),
    (('&Tools',), 'Czjzek &distribution P(C_Q)…', '', False, False),
    (('&Tools',), 'Import &DFT tensors (.magres)…', '', False, False),
    (('&Tools',), 'Re&laxation', '', False, True),
    (('&Tools', 'Re&laxation'), 'Relaxation / series (T1, T2)…', '', False, False),
    (('&Tools', 'Re&laxation'), 'Per-site relaxation…', '', False, False),
    (('&Tools', 'Re&laxation'), 'Variable temperature (Arrhenius / VFT)…', '', False, False),
    (('&Tools',), '&QCPMG', '', False, True),
    (('&Tools', '&QCPMG'), 'QCPMG (echo train → spectrum)…', '', False, False),
    (('&Tools', '&QCPMG'), 'Infinite-field δiso (2 fields)…', '', False, False),
    (('&Tools', '&QCPMG'), 'Batch infinite-field δiso…', '', False, False),
    (('&Tools',), 'R&EDOR (dipolar coupling)…', '', False, False),
    ((), '&View', '', False, True),
    (('&View',), 'Residual', '', True, False),
    (('&View',), 'Components', '', True, False),
    (('&View',), 'Component &labels', '', True, False),
    (('&View',), '&Paddles', '', True, False),
    (('&View',), 'O&verlays', 'Ctrl+Shift+V', True, False),
    (('&View',), 'Clear overlays', '', False, False),
    (('&View',), '&Literature shift ranges', '', True, False),
    (('&View',), 'Scroll &nudges fit values', '', True, False),
    (('&View',), '&Zoom', '', False, True),
    (('&View', '&Zoom'), '&Full spectrum', '', False, False),
    (('&View', '&Zoom'), 'Zoom to &sites', '', False, False),
    (('&View', '&Zoom'), '&Back to 2D map', 'Ctrl+2', False, False),
    (('&View',), 'Axis &unit', '', False, True),
    (('&View', 'Axis &unit'), 'δ (&ppm)', '', True, False),
    (('&View', 'Axis &unit'), '&kHz  (offset from 0 ppm)', '', True, False),
    (('&View', 'Axis &unit'), '&MHz  (offset from 0 ppm)', '', True, False),
    (('&View',), '&Y axis', '', False, True),
    (('&View', '&Y axis'), '&Raw intensity', '', True, False),
    (('&View', '&Y axis'), 'Normalise to &maximum', '', True, False),
    (('&View', '&Y axis'), 'Normalise to &area', '', True, False),
    (('&View', '&Y axis'), 'Normalise to area of a &region…', '', True, False),
    (('&View',), 'Czjzek &width display', '', False, True),
    (('&View', 'Czjzek &width display'), 'σ(Cq)  —  the fitted Czjzek distribution width (mrsimulator σ)', '', True, False),
    (('&View', 'Czjzek &width display'), "Cq ≈ 2σ  —  dmfit's sCZ_CQ = the Czjzek-paper σ_Cz (2σ)", '', True, False),
    (('&View', 'Czjzek &width display'), "CQ (dmfit) = 4σ  —  what dmfit's CQ box displays (2×sCZ_CQ); ≈ the mode of the |Cq| distribution (3.7σ)", '', True, False),
    (('&View', 'Czjzek &width display'), 'P_Q = 2√5·σ  —  rms quadrupolar product √⟨P_Q²⟩ = √5·σ_Cz — field-independent invariant (Edén 2023 Eq. 45)', '', True, False),
    (('&View',), '&Panels', '', False, True),
    (('&View', '&Panels'), 'Explorer', '', True, False),
    (('&View', '&Panels'), 'Datasets', '', True, False),
    (('&View', '&Panels'), 'Workspaces', '', True, False),
    (('&View', '&Panels'), 'Fit parameters', '', True, False),
    (('&View', '&Panels'), 'Report', '', True, False),
    (('&View', '&Panels'), 'Processing', '', True, False),
    (('&View', '&Panels'), 'Fit &health strip', '', True, False),
    (('&View',), '&Theme', '', False, True),
    (('&View', '&Theme'), 'Light', '', True, False),
    (('&View', '&Theme'), 'Sepia (paper)', '', True, False),
    (('&View', '&Theme'), 'Solarized Light', '', True, False),
    (('&View', '&Theme'), 'High Contrast Light', '', True, False),
    (('&View', '&Theme'), 'Dark', '', True, False),
    (('&View', '&Theme'), 'Slate', '', True, False),
    (('&View', '&Theme'), 'Nord', '', True, False),
    (('&View', '&Theme'), 'Ocean', '', True, False),
    (('&View', '&Theme'), 'Solarized Dark', '', True, False),
    (('&View', '&Theme'), 'High Contrast Dark', '', True, False),
    (('&View', '&Theme'), 'More styles…', '', False, True),
    (('&View', '&Theme', 'More styles…'), 'Normal', '', True, False),
    (('&View', '&Theme', 'More styles…'), 'Y2K', '', True, False),
    (('&View', '&Theme', 'More styles…'), 'Dreamcore', '', True, False),
    (('&View', '&Theme', 'More styles…'), 'Gen X Soft Club', '', True, False),
    (('&View', '&Theme', 'More styles…'), 'Vaporwave', '', True, False),
    (('&View',), 'Text &size', '', False, True),
    (('&View', 'Text &size'), 'Small  (8 pt)', '', True, False),
    (('&View', 'Text &size'), 'Normal  (9 pt)', '', True, False),
    (('&View', 'Text &size'), 'Large  (11 pt)', '', True, False),
    (('&View', 'Text &size'), 'Larger  (13 pt)', '', True, False),
    ((), '&Plotting', '', False, True),
    (('&Plotting',), 'Plotting &studio…', '', False, False),
    (('&Plotting',), 'Plot &current spectrum…', '', False, False),
    (('&Plotting',), 'New &2D contour plot…', '', False, False),
    ((), '&Help', '', False, True),
    (('&Help',), '&Command palette…', 'Ctrl+Shift+P', False, False),
    (('&Help',), 'User &manuals', '', False, True),
    (('&Help', 'User &manuals'), 'Getting started', '', False, False),
    (('&Help', 'User &manuals'), '1D spectra — processing & fitting', '', False, False),
    (('&Help', 'User &manuals'), 'Lineshapes — models & physics', '', False, False),
    (('&Help', 'User &manuals'), 'Fitting glasses for publication (Edén 2023)', '', False, False),
    (('&Help', 'User &manuals'), 'Literature shift ranges — data & sources', '', False, False),
    (('&Help', 'User &manuals'), 'Processing reference', '', False, False),
    (('&Help', 'User &manuals'), '2D processing', '', False, False),
    (('&Help', 'User &manuals'), 'MQMAS (2D)', '', False, False),
    (('&Help', 'User &manuals'), 'HMQC & correlation', '', False, False),
    (('&Help', 'User &manuals'), 'Relaxation (T1/T2)', '', False, False),
    (('&Help', 'User &manuals'), 'QCPMG', '', False, False),
    (('&Help', 'User &manuals'), 'Multi-dataset & co-fitting', '', False, False),
    (('&Help', 'User &manuals'), 'DFT tensors — import & shift calibration', '', False, False),
    (('&Help',), '&Tutorials', '', False, True),
    (('&Help', '&Tutorials'), '1 · A first fit — ²⁷Al Czjzek', '', False, False),
    (('&Help', '&Tutorials'), '2 · Constraints — fix, bound, link', '', False, False),
    (('&Help', '&Tutorials'), '3 · Plotting studio figures', '', False, False),
    (('&Help', '&Tutorials'), '4 · Batch fitting a composition series', '', False, False),
    (('&Help', '&Tutorials'), '5 · MQMAS — 2D processing & fitting', '', False, False),
    (('&Help', '&Tutorials'), '6 · Error analysis — covariance, Monte-Carlo, χ² profile', '', False, False),
    (('&Help', '&Tutorials'), '7 · Static wideline — ⁸¹Br WURST-CPMG', '', False, False),
    (('&Help',), 'About LARMOR', '', False, False),
    (('&Help',), 'More…', '', False, False),
)

SHORTCUTS = frozenset(['Ctrl+2', 'Ctrl+Alt+S', 'Ctrl+F', 'Ctrl+I', 'Ctrl+O', 'Ctrl+P', 'Ctrl+S', 'Ctrl+Shift+A', 'Ctrl+Shift+C', 'Ctrl+Shift+D', 'Ctrl+Shift+E', 'Ctrl+Shift+O', 'Ctrl+Shift+P', 'Ctrl+Shift+S', 'Ctrl+Shift+V', 'Ctrl+T', 'Ctrl+Y', 'Ctrl+Z', 'F5', 'F6', 'F7', 'F9'])
# main toolbar: Undo, Redo, the '＋ Add line' split button (a widget action),
# five quick model buttons, the placing label (widget action) = 9; the
# sidebar's 9 view actions (larmor/desktop/addline_toolbar.py, wip/WB)
TOOLBAR_ACTIONS = 18
# QMenu count under the menu bar: 35


# ------------------------------------------------------------------ helpers
def _mixin_files():
    return sorted(DESKTOP.glob("mw_*.py"))


def _mixins():
    from PySide6.QtWidgets import QMainWindow

    from larmor.desktop.app import MainWindow

    return [b for b in MainWindow.__bases__ if b is not QMainWindow]


def _defined(cls):
    """Names a class body defines itself: methods, staticmethods, constants;
    ``__init__`` kept, other dunders and PySide's staticMetaObject dropped."""
    return {k for k in vars(cls)
            if (not k.startswith("__") or k == "__init__")
            and k != "staticMetaObject"}


def _menu_rows(win):
    rows = []

    def walk(menu, path):
        for a in menu.actions():
            if a.isSeparator():
                continue
            sub = a.menu()
            rows.append((path, a.text(), a.shortcut().toString(),
                         a.isCheckable(), sub is not None))
            if sub is not None and a.text() not in EXCLUDE_CHILDREN:
                walk(sub, path + (a.text(),))

    for a in win.menuBar().actions():
        sub = a.menu()
        rows.append(((), a.text(), a.shortcut().toString(), a.isCheckable(),
                     sub is not None))
        if sub is not None:
            walk(sub, (a.text(),))
    return rows


@pytest.fixture()
def win(monkeypatch):
    from PySide6.QtWidgets import QApplication

    monkeypatch.setenv("LARMOR_NO_SESSION", "1")   # never touch a real session
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    yield w
    w.close()


# ------------------------------------------------------- static contract
def test_public_surface_unchanged():
    """Every pre-split member is still defined on MainWindow or on one of its
    mixins -- never merely inherited from QMainWindow under the same name."""
    from larmor.desktop.app import MainWindow

    assert len(SURFACE) == len(set(SURFACE))
    owners = [MainWindow] + _mixins()
    missing = [n for n in SURFACE if not any(n in vars(c) for c in owners)]
    assert not missing, missing


def test_facade_holds_only_construction_and_qt_overrides():
    from larmor.desktop.app import MainWindow

    assert _defined(MainWindow) == FACADE_METHODS


def test_mixins_are_plain_and_disjoint():
    from PySide6.QtCore import QObject, Signal

    mixins = _mixins()
    seen = {}
    for b in mixins:
        assert not issubclass(b, QObject), b
        assert "__init__" not in vars(b), b
        assert not any(isinstance(v, Signal) for v in vars(b).values()), b
        assert not (QT_OVERRIDES & set(vars(b))), (b, QT_OVERRIDES & set(vars(b)))
        for name in _defined(b):
            assert name not in seen, (
                f"{name} is defined on both {seen[name].__name__} and "
                f"{b.__name__}; MRO order would hide one of them")
            seen[name] = b


def test_self_method_calls_resolve():
    """Every ``self.X(`` inside MainWindow or a mixin names a method of the
    assembled class or an attribute assigned on self somewhere."""
    from larmor.desktop.app import MainWindow

    calls, stores = [], set()
    for f in [DESKTOP / "app.py"] + _mixin_files():
        tree = ast.parse(f.read_text(encoding="utf-8"))
        for cls in tree.body:
            if not isinstance(cls, ast.ClassDef):
                continue
            if cls.name != "MainWindow" and not cls.name.endswith("Mixin"):
                continue
            for node in ast.walk(cls):
                if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) \
                        and node.value.id == "self" and isinstance(node.ctx, ast.Store):
                    stores.add(node.attr)
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                        and isinstance(node.func.value, ast.Name) \
                        and node.func.value.id == "self":
                    calls.append((f.name, node.lineno, node.func.attr))
    assert len(calls) > 500, len(calls)          # the scan found the code
    unresolved = [f"{fn}:{ln}  self.{name}(" for fn, ln, name in calls
                  if not hasattr(MainWindow, name) and name not in stores]
    assert not unresolved, "\n".join(unresolved)


def test_mixins_do_not_import_app():
    for f in _mixin_files() + [DESKTOP / "workers.py"]:
        if not f.exists():
            continue
        tree = ast.parse(f.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert node.module != "larmor.desktop.app", (f.name, node.lineno)
                if node.module == "larmor.desktop":
                    assert not any(a.name == "app" for a in node.names), (f.name, node.lineno)
            elif isinstance(node, ast.Import):
                assert not any(a.name == "larmor.desktop.app" for a in node.names), \
                    (f.name, node.lineno)


def test_reexports_still_import():
    import importlib

    app = importlib.import_module("larmor.desktop.app")
    for name in REEXPORTS:
        assert hasattr(app, name), name
    assert set(REEXPORTS) <= set(app.__all__), set(REEXPORTS) - set(app.__all__)
    for name in app.__all__:
        assert hasattr(app, name), name
    for name, home in HOMES.items():
        mod = importlib.import_module(home)
        assert getattr(app, name) is getattr(mod, name), (name, home)


# ------------------------------------------------------- window contract
def test_menu_tree_and_shortcuts_golden(win):
    rows = _menu_rows(win)
    if tuple(rows) != GOLDEN_MENU:
        fresh = "\n".join(f"    {r!r}," for r in rows)
        pytest.fail("the menu tree changed; if that is intended, replace "
                    "GOLDEN_MENU with:\nGOLDEN_MENU = (\n" + fresh + "\n)")

    from PySide6.QtGui import QAction
    from PySide6.QtWidgets import QToolBar

    shortcuts = [a.shortcut().toString() for a in win.findChildren(QAction)
                 if a.shortcut().toString()]
    assert frozenset(shortcuts) == SHORTCUTS, sorted(set(shortcuts) ^ SHORTCUTS)
    dupes = sorted({s for s in shortcuts if shortcuts.count(s) > 1})
    assert not dupes, dupes                      # Qt fires neither on a clash
    n_toolbar = sum(1 for tb in win.findChildren(QToolBar)
                    for a in tb.actions() if not a.isSeparator())
    assert n_toolbar == TOOLBAR_ACTIONS

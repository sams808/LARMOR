"""XFact — a small live-data 'random fact' popup card, organized as
independent packs sharing one Qt dialog, one background-thread worker,
and one HTTP layer. Vendored into LARMOR from github.com/sams808/XFact.

    from xfact import show_fact
    show_fact(parent=my_main_window)             # random pack each time
    show_fact(parent=my_main_window, pack="birds") # pin to one pack

See xfact/core/registry.py to add a new pack.
"""

from .core.dialog import FactDialog
from .core.registry import get_pack, get_random_card_from_any_pack, list_packs


def show_fact(parent=None, pack: str | None = None) -> FactDialog:
    get_card = get_random_card_from_any_pack if pack is None else get_pack(pack).get_random_card
    dlg = FactDialog(get_card, parent)
    dlg.show()
    return dlg


__all__ = ["show_fact", "list_packs"]

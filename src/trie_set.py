"""Deterministic membership for canonical HNP integer/point encodings.

Integers, None, and nested tuples/lists are encoded as compact JSON bytes;
tuple and list encodings of the same point deliberately coincide. Trie
membership and insertion take O(serialized-key length) operations. Integer
serialization itself uses polynomial bit time. Child dictionary keys are
bytes, so every dictionary has at most 256 entries: no unbounded hash-table
collision assumption enters the per-key bound.
"""
from __future__ import annotations

import json


def canonical_bytes(value) -> bytes:
    """Encode the supported public canonical values without ambiguity."""
    def normalized(item):
        if item is None:
            return None
        if type(item) is int:
            return item
        if isinstance(item, (tuple, list)):
            return [normalized(part) for part in item]
        raise TypeError("TrieSet keys must be integers, None, tuples, or lists")

    return json.dumps(normalized(value), separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False).encode("ascii")


class _Node:
    __slots__ = ("children", "terminal")

    def __init__(self):
        self.children: dict[int, _Node] = {}
        self.terminal = False


class TrieSet:
    """Small set interface with deterministic work per canonical encoding."""
    __slots__ = ("_root", "_size")

    def __init__(self, values=()):
        self._root = _Node()
        self._size = 0
        for value in values:
            self.add(value)

    def __len__(self):
        return self._size

    def __contains__(self, value):
        node = self._root
        for octet in canonical_bytes(value):
            node = node.children.get(octet)
            if node is None:
                return False
        return node.terminal

    def add(self, value) -> bool:
        """Insert a key; return True exactly when its encoding was absent."""
        node = self._root
        for octet in canonical_bytes(value):
            child = node.children.get(octet)
            if child is None:
                child = _Node()
                node.children[octet] = child
            node = child
        if node.terminal:
            return False
        node.terminal = True
        self._size += 1
        return True

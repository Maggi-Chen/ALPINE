"""
Shared module for column ordering in variant count tables.

This module defines the fixed column order for variant classification tables
and provides utility functions for sorting columns consistently across
both per-sample counting and multi-sample merging operations.
"""

# Fixed column order: exact matches and prefixes for transgene-suffixed categories
# Prefixes end with '-' to indicate they match categories like 'HDR-CD19', 'HDR-HLAE'
COLUMN_ORDER = [
    'Unmodified',
    'Unmodified-with-SNP',
    'DEL-small',
    'INS-small',
    'DEL-large',
    'INS-large',
    'HDR-',                      # prefix: HDR-{Transgene}
    'Non-HDR-with-ITR-',         # prefix: Non-HDR-with-ITR-{Transgene}
    'Non-HDR-without-ITR-',      # prefix: Non-HDR-without-ITR-{Transgene}
    'INV',
    'DUP',
    'Unclassified',
]


def get_sort_key(column_name):
    """
    Return a sort key tuple (order_index, suffix) for a column name.

    Args:
        column_name: The column name to generate a sort key for.

    Returns:
        tuple: (order_index, suffix) where:
            - order_index: position in COLUMN_ORDER (lower = earlier)
            - suffix: for prefix matches, enables alphabetical sorting within group
    """
    for idx, pattern in enumerate(COLUMN_ORDER):
        if pattern.endswith('-'):
            # Prefix match (e.g., 'HDR-' matches 'HDR-CD19')
            if column_name.startswith(pattern):
                suffix = column_name[len(pattern):]
                return (idx, suffix)
        else:
            # Exact match
            if column_name == pattern:
                return (idx, '')
    # Unknown categories go to the end, sorted alphabetically
    return (len(COLUMN_ORDER), column_name)


def sort_columns_with_counts(read_count):
    """
    Sort columns according to COLUMN_ORDER and add missing exact-match categories with count 0.

    This function is used by count_classification_table.py for per-sample counting.

    Args:
        read_count: Dictionary mapping column names to counts.

    Returns:
        list: Ordered list of (column_name, count) tuples.
    """
    existing_columns = set(read_count.keys())
    ordered_columns = []

    for pattern in COLUMN_ORDER:
        if pattern.endswith('-'):
            # Prefix pattern: find all matching columns, sort alphabetically by suffix
            matching = [col for col in existing_columns if col.startswith(pattern)]
            matching.sort(key=lambda x: x[len(pattern):])  # sort by suffix
            for col in matching:
                ordered_columns.append((col, read_count[col]))
                existing_columns.discard(col)
        else:
            # Exact match pattern
            if pattern in existing_columns:
                ordered_columns.append((pattern, read_count[pattern]))
                existing_columns.discard(pattern)
            else:
                # Missing category: add with count 0
                ordered_columns.append((pattern, 0))

    # Append any remaining unknown columns at the end (sorted alphabetically)
    for col in sorted(existing_columns):
        ordered_columns.append((col, read_count[col]))

    return ordered_columns


def sort_columns(all_tags):
    """
    Sort column names according to COLUMN_ORDER and add missing exact-match categories.

    This function is used by merge_count_table.py for merging multiple samples.

    Args:
        all_tags: Collection of column names from all samples.

    Returns:
        list: Ordered list of column names.
    """
    existing_tags = set(all_tags)
    ordered_columns = []

    for pattern in COLUMN_ORDER:
        if pattern.endswith('-'):
            # Prefix pattern: find all matching columns, sort alphabetically by suffix
            matching = [tag for tag in existing_tags if tag.startswith(pattern)]
            matching.sort(key=lambda x: x[len(pattern):])  # sort by suffix
            for tag in matching:
                ordered_columns.append(tag)
                existing_tags.discard(tag)
        else:
            # Exact match pattern: always include
            ordered_columns.append(pattern)
            existing_tags.discard(pattern)

    # Append any remaining unknown columns at the end (sorted alphabetically)
    for tag in sorted(existing_tags):
        ordered_columns.append(tag)

    return ordered_columns

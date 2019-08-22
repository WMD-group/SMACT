"""Tools for database interfacing for high throughput IO."""

from typing import Sequence, List

from .structure import SmactStructure
import sqlite3


class StructureDB:
    """SQLite Structure Database interface.

    Acts as a context manager for database interfacing
    and wraps several useful SQLite commands within
    methods.

    Attributes:
        db: The database name.
        conn: The database connection. Only open when
            used as a context manager.
        cur: The database connection cursor. Only usable
            when class implemented as context manager.

    Examples:
        Connecting to a database in memory:
        >>> DB = StructureDB(':memory:')
        >>> with DB as c:
        ...     _ = c.execute("CREATE TABLE test (id, val)")
        ...     c.execute("SELECT * FROM test").fetchall()
        []
        >>> DB.cur.execute("SELECT * FROM test").fetchall()
        Traceback (most recent call last):
            ...
        sqlite3.ProgrammingError: Cannot operate on a closed database.

    """

    def __init__(self, db: str):
        """Set database name."""
        self.db = db

    def __enter__(self) -> sqlite3.Cursor:
        """Initialize database connection."""
        self.conn = sqlite3.connect(self.db)
        self.cur = self.conn.cursor()

        return self.cur

    def __exit__(self, *args):
        """Close database connection."""
        self.conn.commit()
        self.conn.close()

    def add_table(self, table: str) -> bool:
        """Add a table to the database.

        Args:
            table: The name of the table to add

        Returns:
            bool: Whether the operation was successful.

        """
        with self as c:
            c.execute(
              f"""CREATE TABLE {table}
                (composition TEXT NOT NULL, structure TEXT NOT NULL)""",
            )
        return True

    def add_struct(self, struct: SmactStructure, table: str) -> bool:
        """Add a SmactStructure to a table.

        Args:
            struct: The :class:`~.SmactStructure` to add.
            table: The name of the table to add the structure to.

        Returns:
            bool: Whether the operation was successful.

        """
        entry = (struct.composition(), struct.as_poscar())

        with self as c:
            c.execute(f"INSERT into {table} VALUES (?, ?)", entry)

        return True

    def add_structs(self, structs: Sequence[SmactStructure], table: str) -> bool:
        """Add several SmactStructures to a table.

        Args:
            structs: Iterable of :class:`~.SmactStructure`s to add to table.
            table: The name of the table to add the structs to.

        Returns:
            bool: Whether the operation was successful.

        """
        with self as c:
            entries = [(
              struct.composition(),
              struct.as_poscar(), ) for struct in structs]

            c.executemany(f"INSERT into {table} VALUES (?, ?)", entries)

        return True

    def get_structs(self, composition: str, table: str) -> List[SmactStructure]:
        """Get SmactStructures for a given composition.

        Args:
            composition: The composition to search for.
                See :meth:`SmactStructure.composition`.
            table: The name of the table in which to search.

        Returns:
            A list of :class:`~.SmactStructure` s.

        """
        with self as c:
            c.execute(
              f"SELECT structure FROM {table} WHERE composition = ?",
              (composition, ), )
            structs = c.fetchall()
        return [SmactStructure.from_poscar(pos[0]) for pos in structs]

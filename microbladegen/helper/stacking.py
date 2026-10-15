import enum

class StackType(enum.Enum):
    """class defining the type of stacking for Airfoil2D profiles

    Args:
        enum (enum.Emum): inherits enum
    """
    Leading_Edge = 1
    Centroid = 2
    Trailing_Edge = 3
    None_ = 4
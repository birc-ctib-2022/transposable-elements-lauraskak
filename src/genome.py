"""A circular genome for simulating transposable elements."""

from abc import (
    # A tag that says that we can't use this class except by specialising it
    ABC,
    # A tag that says that this method must be implemented by a child class
    abstractmethod
)


class Genome(ABC):
    """Representation of a circular Genome."""

    def __init__(self, n: int):
        """Create a genome of size n."""
        ...  # not implemented yet

    @abstractmethod
    def insert_te(self, pos: int, length: int) -> int:
        """
        Insert a new transposable element.

        Insert a new transposable element at position pos and len
        nucleotide forward.

        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.

        Returns a new ID for the transposable element.
        """
        ...  # not implemented yet

    @abstractmethod
    def copy_te(self, te: int, offset: int) -> int | None:
        """
        Copy a transposable element.

        Copy the transposable element te to an offset from its current
        location.

        The offset can be positive or negative; if positive the te is copied
        upwards and if negative it is copied downwards. If the offset moves
        the copy left of index 0 or right of the largest index, it should
        wrap around, since the genome is circular.

        If te is not active, return None (and do not copy it).
        """
        ...  # not implemented yet

    @abstractmethod
    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        ...  # not implemented yet

    @abstractmethod
    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        ...  # not implemented yet

    @abstractmethod
    def __len__(self) -> int:
        """Get the current length of the genome."""
        ...  # not implemented yet

    @abstractmethod
    def __str__(self) -> str:
        """
        Return a string representation of the genome.

        Create a string that represents the genome. By nature, it will be
        linear, but imagine that the last character is immidiatetly followed
        by the first.

        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        ...  # not implemented yet


class ListGenome(Genome):
    """
    Representation of a genome.

    Implements the Genome interface using Python's built-in lists
    """

    def __init__(self, n: int): #DONE
        """Create a new genome with length n."""
        self.genome_list = ["-"] * n
        self.TE_dict = {}
            
    def insert_te(self, pos: int, TE_length: int) -> int: #DONE
        """
        Insert a new transposable element.

        Insert a new transposable element at position pos and len
        nucleotide forward.

        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.

        Returns a new ID for the transposable element.
        """
        TE = ["A"] * TE_length
        
        # The following checks if the TE insertion collides with an existing active TE
        if self.genome_list[pos] == "A" and self.genome_list[pos+1] == "A":
            for TE_id in self.TE_dict:
                #The following identifies the active TE to be disabled and disables it.
                if self.TE_dict[TE_id][0] < pos and pos < self.TE_dict[TE_id][1]:
                    self = self.disable_te(TE_id) 
                    break
                    
        # the following updates the genome list to include the TE
        self.genome_list = self.genome_list[:pos] + TE + self.genome_list[pos:]
        
        #updates TE_dict
        for TE_id in self.TE_dict:
            if self.TE_dict[TE_id][0] > pos:
                self.TE_dict[TE_id][0] += TE_length
                self.TE_dict[TE_id][1] += TE_length
        
        ## Inputs the new TE in the TE_id dict with the start and end positions. 
        previous_max_id = max(self.TE_dict.keys)
        self.TE_dict[previous_max_id+1] = [pos + 1, pos + TE_length]
        
        return previous_max_id+1

    def copy_te(self, TE_id: int, offset: int) -> int | None: #DONE
        """
        Copy a transposable element.

        Copy the transposable element te to an offset from its current
        location.

        The offset can be positive or negative; if positive the te is copied
        upwards and if negative it is copied downwards. If the offset moves
        the copy left of index 0 or right of the largest index, it should
        wrap around, since the genome is circular.

        If te is not active, return None (and do not copy it).
        """
        
        if self.TE_dict[TE_id] == None:
            return
        else:
            start = self.TE_dict[TE_id][0]
            end = self.TE_dict[TE_id][1]
            TE_length = end - start + 1
            genome_length = len(self)
            
            if offset > 0:
                ## Situation where the offset is positive
                pos = (end + offset) % genome_length
                self.insert_te(pos, TE_length)
            
            elif offset < 0:
                ## Situation where the offset is negative
                pos = genome_length + start + offset - 2
                
            else: 
                pos = end
                self.insert_te(pos, TE_length)
                
    def disable_te(self, TE_id: int) -> None: #DONE
        """ Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        #lav A'er om til X'er
        start = self.TE_dict[TE_id][0]
        end = self.TE_dict[TE_id][1]
        TE_length = end - start + 1
        
        inactive_TE = ["X"] * TE_length
        
        # Updates the genome list
        self.genome_list = self.genome_list[:start] + inactive_TE + self.genome_list[end+1:]
        
        # Inactivates the TE in the TE_id dict
        self.TE_dict[TE_id] = None

    def active_tes(self) -> list[int]: # DONE
        """Get the active TE IDs."""
        
        active_TEs = []
        
        for TE_id in self.TE_dict:
            if self.TE_dict[TE_id] != None:
                active_TEs.append(TE_id)
                
        return active_TEs

    def __len__(self) -> int: #DONE
        """Current length of the genome."""
        return len(self.genome_list)

    def __str__(self) -> str: #DONE
        """
        Return a string representation of the genome.

        Create a string that represents the genome. By nature, it will be
        linear, but imagine that the last character is immidiatetly followed
        by the first.

        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        return "".join(self.genome_list)


class Link_element():
    def __init__(self, nucleotide, next):
        self.nucleotide = nucleotide
        self.next = next
        
            

class LinkedListGenome(Genome):
    """
    Representation of a genome.

    Implements the Genome interface using linked lists.
    """
    


    def __init__(self, n: int):
        """Create a new genome with length n."""
        
        self.head = Link_element("-", None)
        
        next_elem = self.head
        count = n-1
        
        while count > 0:
            link = Link_element("-", next_elem)
            
            next_elem = link
            count -= 1
        
        (self.head).next = next_elem
        

    def insert_te(self, pos: int, length: int) -> int:
        """
        Insert a new transposable element.

        Insert a new transposable element at position pos and len
        nucleotide forward.

        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.

        Returns a new ID for the transposable element.
        """
        
         
        ...  # FIXME
        return -1

    def copy_te(self, te: int, offset: int) -> int | None:
        """
        Copy a transposable element.

        Copy the transposable element te to an offset from its current
        location.

        The offset can be positive or negative; if positive the te is copied
        upwards and if negative it is copied downwards. If the offset moves
        the copy left of index 0 or right of the largest index, it should
        wrap around, since the genome is circular.

        If te is not active, return None (and do not copy it).
        """
        ...  # FIXME

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        ...  # FIXME

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        # FIXME
        return []

    def __len__(self) -> int:
        """Current length of the genome."""
        # FIXME
        return 0

    def __str__(self) -> str:
        """
        Return a string representation of the genome.

        Create a string that represents the genome. By nature, it will be
        linear, but imagine that the last character is immidiatetly followed
        by the first.

        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        return "FIXME"


DNA = ListGenome(10)

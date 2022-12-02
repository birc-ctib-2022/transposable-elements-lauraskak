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
        
        return self.insert_te(pos, length)
    

    @abstractmethod
    def copy_te(self, te: int, offset: int):
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
        
        self.copy_te(te, offset)
        

    @abstractmethod
    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        
        self.disable_te(te)
        

    @abstractmethod
    def active_tes(self) -> list:
        """Get the active TE IDs."""
        
        return self.active_tes()
        

    @abstractmethod
    def __len__(self) -> int:
        """Get the current length of the genome."""
        return len(self)

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
        return str(self)


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
        assert pos >= 1
        
        TE = ["A"] * TE_length
        
        genome_length = len(self)
        
        # Makes sure that pos is a position in the genome.
        pos = pos % genome_length - 1
        
        # The following checks if the TE insertion collides with an existing active TE
        if self.genome_list[pos] == "A" and self.genome_list[pos+1] == "A":
            for TE_id in self.TE_dict:
                #The following identifies the active TE to be disabled and disables it.
                
                if self.TE_dict[TE_id] != None:
                    if self.TE_dict[TE_id][0] <= pos and pos < self.TE_dict[TE_id][1] and self.TE_dict[TE_id] != None:
                        self.disable_te(TE_id)
                        break
                
        # The following updates the genome list to include the TE
        self.genome_list = self.genome_list[:pos+1] + TE + self.genome_list[pos+1:]
        
        
        # Updates TE_dict
        for TE_id in self.TE_dict:
            if self.TE_dict[TE_id] != None:
                if self.TE_dict[TE_id][0] > pos:
                    self.TE_dict[TE_id][0] += TE_length
                    self.TE_dict[TE_id][1] += TE_length
    
        
        ## Inputs the new TE in the TE_id dict with the start and end positions. 
        if self.TE_dict == {}:
            new_id = 1
        else:
            previous_max_id = max(self.TE_dict.keys())
            new_id = previous_max_id + 1
            
        self.TE_dict[new_id] = [pos + 1, pos + TE_length]
        
        return new_id
    

    def copy_te(self, TE_id: int, offset: int): #DONE
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
                ## Situation where the offset is positive, so the TE should be moved downstream
                
                if start + offset >= genome_length:
                    pos = (start + offset) % genome_length - 1
                    new_id = self.insert_te(pos, TE_length)
                else: 
                    pos = start + offset
                    new_id = self.insert_te(pos, TE_length)
            
            elif offset < 0:
                ## Situation where the offset is negative, so the TE should be moved upstream
                
                if start-1 <= -offset:
                    diff = - offset - start
                    pos = genome_length - diff
                    new_id = self.insert_te(pos, TE_length)
                else:
                    pos = start + offset - 1
                    new_id = self.insert_te(pos, TE_length)
                
            else: 
                # If the offset is 0
                pos = end
                new_id = self.insert_te(pos, TE_length)
                
            return new_id
                
    def disable_te(self, TE_id: int) -> None: #DONE
        """ Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        if TE_id not in self.TE_dict:
            # Makes sure it does nothing if the TE_id doesn't exist
            return
        
        else:
            #lav A'er om til X'er
            start = self.TE_dict[TE_id][0]
            end = self.TE_dict[TE_id][1]
            TE_length = end - start + 1
            
            inactive_TE = ["x"] * TE_length
            
            # Updates the genome list
            self.genome_list = self.genome_list[:start] + inactive_TE + self.genome_list[end+1:]
            
            # Inactivates the TE in the TE_id dict
            self.TE_dict[TE_id] = None

    def active_tes(self) -> list: # DONE
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
    def __init__(self, previous, next):
        self.next = next
        self.previous = previous
        self.is_TE = False


class LinkedListGenome(Genome):
    """
    Representation of a genome.

    Implements the Genome interface using linked lists.
    """
    


    def __init__(self, n: int):
        """Create a new genome with length n."""
        assert n > 0       
        
        self.head = Link_element(None, None)
        
        previous_elem = self.head
        link = self.head
        
        count = n-1
        
        while count > 0:
            link = Link_element(previous_elem, None)
            previous_elem.next = link
            previous_elem = link
            count -= 1
        
        self.head.previous = link
        link.next = self.head
        

    def insert_te(self, pos: int, TE_length: int) -> int:
        """
        Insert a new transposable element.

        Insert a new transposable element at position pos and len
        nucleotide forward.

        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.

        Returns a new ID for the transposable element.
        """
        
        # What should the new TE_id be?
        
        link = self.head
        
        max_id = 0
        
        looped = False
        
        while not looped:
            
            if link.is_TE: 
                if link.TE_id > max_id:
                    max_id = link.TE_id
            
            if link.next == self.head:
                looped = True
                    
            link = link.next
        
        new_id = max_id + 1
        
        # Now max_id is the maximal id count and new_id is the name for the TE to be inserted
        
        link = self.head
        
        pos_count = 1
        
        while pos_count != pos:
      
            link = link.next
            pos_count += 1 
            
        # Now the link is the position which the TE should be inserted after
        
        
        # Disable if the TE is about to be put into the middle of another TE
        
        if link.is_TE and link.next.is_TE:
            if link.TE_id == link.next.TE_id:
                self.disable_te(link.TE_id)
        
        ## Inserts the new TE
        previous_elem = link
        old_next = link.next
        
        TE_count = 0
        
        while TE_count != TE_length:
            
            link = Link_element(previous_elem, None)
            link.is_TE = True
            link.TE_id = new_id
            link.is_active = True
            previous_elem.next = link
            previous_elem = link
            
            TE_count += 1
        
        link.next = old_next
        
        return new_id
        
        

    def copy_te(self, TE_id: int, offset: int):
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
        
        ## Finds out how long the TE is
        
        link = self.head
        
        looped = False
        
        TE_length = 0
        
        pos_count = 0
        
        TE_start_pos = None
        
        genome_length = len(self)
                
        while not looped:
            if link.is_TE:
                if link.TE_id == TE_id:
                    if link.is_active:
                        if TE_start_pos == None:
                            # Marks the TE's start position
                            TE_start_pos = pos_count
                        TE_length += 1
                    else:
                        # Stops the operation if the TE about to be copied is inaktive
                        return
            
            if link.next == self.head:
                looped = True
            
            link = link.next
            
            pos_count += 1
            
        
        # now TE_length should be the length of the TE to be copied.
        # and TE start and end is defined.
        
        ## If the offset is positive
        
        if offset >= 0:
            insert_pos = TE_start_pos + offset % genome_length
            new_id = self.insert_te(insert_pos, TE_length)
            
        ## If the offset is negative or 0
        
        else:
            if -offset >= TE_start_pos - 1:
                diff = - offset - TE_start_pos
                insert_pos = genome_length - diff
                new_id = self.insert_te(insert_pos, TE_length)
            else:
                insert_pos = TE_start_pos - offset - 1
                print(insert_pos, TE_length)
                new_id = self.insert_te(insert_pos, TE_length)
                
        
        return new_id
            

    def disable_te(self, TE_id: int): #DONE
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        link = self.head
        
        looped = False
        
        while not looped:
            
            if link.is_TE:
                if link.TE_id == TE_id:
                    link.is_active = False
            
            if link.next == self.head:
                looped = True
                
            link = link.next       
    
    def active_tes(self) -> list : #DONE
        """Get the active TE IDs."""
        
        link = self.head
        ## Since self.head cannot be a TE, we can ignore this
        
        last_TE_id = None
        active_TE_list = []
        
        looped = False
        
        while not looped:
            
            if link.is_TE:
                if link.TE_id != last_TE_id:
                    
                    if link.is_active:
                        active_TE_list.append(link.TE_id)
                        
                    last_TE_id = link.TE_id
            
            if link.next == self.head:
                looped = True
            
            link = link.next
                
        return active_TE_list

    def __len__(self) -> int: #DONE
        """Current length of the genome."""
        
        link = self.head
        count = 1
        
        while link.next != self.head:
            link = link.next
            count += 1
            
        return count

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
        
        genome = ""
        
        link = self.head
        
        looped = False
        
        while not looped:
            
            if link.is_TE:
                if link.is_active:
                    genome += "A"
                else:
                    genome += "x"
            else:
                genome += "-"
                
            if link.next == self.head:
                looped = True
            
            link = link.next
            
        return genome


### CHECK OF ListGenome functions

# DNA = ListGenome(10)

# print("Genome list: ", DNA.genome_list)
# print("TE dict: ", DNA.TE_dict)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("inserts 3 in after pos 4")
# DNA.insert_te(2, 3)

# print("inserts 5 in after pos 4")
# DNA.insert_te(4, 5)

# print("inserts copy 5 after TE 2")
# DNA.copy_te(2, 5)

# print("Genome list: ", DNA.genome_list)
# print("TE dict: ", DNA.TE_dict)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("expected: ")
# print("---XXAAAAAX----AAAAA---")
# print("result: ")
# print(str(DNA))


### CHECK OF LinkedListGenome functions

# DNA = LinkedListGenome(10)

# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("inserts 3 after pos 6")
# DNA.insert_te(6, 3)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("Disable TE id that does not exist")
# DNA.disable_te(3)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("inserts 4 after pos 7")
# DNA.insert_te(7, 4)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("Prints the active TE ids")
# print(DNA.active_tes())

# print("inserts TE of length 6 after pos 15, so achually at pos 21 or pos 4")
# DNA.insert_te(21, 6)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("copies id 2 with offset of 5")
# DNA.copy_te(2, 5)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("Prints the active TE ids")
# print(DNA.active_tes())


# print("copies id 4 with offset of -2")
# DNA.copy_te(4, -2)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))

# print("inserts TE of length 2 at pos 6")
# DNA.insert_te(6, 2)
# print("Length of the seq: ", len(DNA))
# print("print of the sequence: ", str(DNA))


### Thomas' test

genome = LinkedListGenome(20)
#genome = ListGenome(20)

assert str(genome) == "--------------------"
assert genome.active_tes() == []

assert 1 == genome.insert_te(5, 10)   # Insert te 1
assert str(genome) == "-----AAAAAAAAAA---------------"
assert genome.active_tes() == [1]

assert 2 == genome.insert_te(10, 10)  # Disable 1 but make 2 active
assert str(genome) == "-----xxxxxAAAAAAAAAAxxxxx---------------"
assert genome.active_tes() == [2]

# Make TE 3 20 to the right of the start of 2
assert 3 == genome.copy_te(2, 20)
assert str(genome) == "-----xxxxxAAAAAAAAAAxxxxx-----AAAAAAAAAA----------"
assert genome.active_tes() == [2, 3]

# Make TE 4 15 to the leftt of the start of 2
assert 4 == genome.copy_te(2, -15)
assert str(genome) == "-----xxxxxAAAAAAAAAAxxxxx-----AAAAAAAAAA-----AAAAAAAAAA-----"
assert genome.active_tes() == [2, 3, 4]

assert 5 == genome.insert_te(50, 10)
assert str(genome) == "-----xxxxxAAAAAAAAAAxxxxx-----AAAAAAAAAA-----xxxxxAAAAAAAAAAxxxxx-----"
assert genome.active_tes() == [2, 3, 5]

genome.disable_te(3)
assert str(genome) == "-----xxxxxAAAAAAAAAAxxxxx-----xxxxxxxxxx-----xxxxxAAAAAAAAAAxxxxx-----"
assert genome.active_tes() == [2, 5]






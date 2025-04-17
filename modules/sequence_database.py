"""
Sequence Database module for fetching and managing DNA sequences from various repositories
"""
import os
import logging
import requests
import io
import time
import re
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from typing import Optional, Dict, List, Union, Tuple
from io import StringIO

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configure Entrez
Entrez.email = "your_email@example.com"  # Replace with a project email

class SequenceRepository:
    """Base class for sequence repositories"""
    
    def fetch_sequence(self, identifier: str) -> Optional[SeqRecord]:
        """Fetch a sequence from the repository"""
        raise NotImplementedError("Subclasses must implement this method")
    
    def search_sequence(self, query: str, max_results: int = 10) -> List[Dict]:
        """Search for sequences matching the query"""
        raise NotImplementedError("Subclasses must implement this method")

class NCBIRepository(SequenceRepository):
    """Repository for fetching sequences from NCBI"""
    
    def fetch_sequence(self, identifier: str) -> Optional[SeqRecord]:
        """
        Fetch a sequence from NCBI by accession number or GI
        
        Args:
            identifier: NCBI accession number or GI
            
        Returns:
            SeqRecord object or None if not found
        """
        try:
            logger.info(f"Fetching sequence from NCBI: {identifier}")
            
            # Handle rate limiting (max 3 requests per second)
            time.sleep(0.33)
            
            # Try to fetch the sequence
            handle = Entrez.efetch(db="nucleotide", id=identifier, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            # Ensure molecule_type is set
            if not hasattr(record, 'annotations') or not record.annotations:
                record.annotations = {}
            if 'molecule_type' not in record.annotations:
                record.annotations['molecule_type'] = 'DNA'
                
            logger.info(f"Successfully fetched {identifier} from NCBI ({len(record.seq)} bp)")
            return record
            
        except Exception as e:
            logger.error(f"Failed to fetch {identifier} from NCBI: {str(e)}")
            return None
    
    def search_sequence(self, query: str, max_results: int = 10) -> List[Dict]:
        """
        Search for sequences in NCBI matching the query
        
        Args:
            query: Search term
            max_results: Maximum number of results to return
            
        Returns:
            List of dictionaries with sequence metadata
        """
        try:
            logger.info(f"Searching NCBI for: {query}")
            
            # Handle rate limiting
            time.sleep(0.33)
            
            # Perform the search
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()
            
            # If no results, return empty list
            if not record["IdList"]:
                logger.info(f"No results found for query: {query}")
                return []
            
            # Fetch summary for each result
            id_list = ",".join(record["IdList"])
            handle = Entrez.esummary(db="nucleotide", id=id_list)
            summaries = Entrez.read(handle)
            handle.close()
            
            # Format results
            results = []
            for summary in summaries:
                results.append({
                    "id": summary["AccessionVersion"],
                    "title": summary.get("Title", ""),
                    "length": summary.get("Length", 0),
                    "source": "NCBI",
                    "type": "DNA",
                    "description": f"{summary.get('Title', '')} [{summary.get('Organism', '')}]"
                })
            
            logger.info(f"Found {len(results)} results for query: {query}")
            return results
            
        except Exception as e:
            logger.error(f"Error searching NCBI: {str(e)}")
            return []

class AddgeneRepository(SequenceRepository):
    """Repository for fetching sequences from Addgene"""
    
    def __init__(self):
        self.base_url = "https://www.addgene.org/search/sequence/"
        self.api_url = "https://www.addgene.org/api/v1"
    
    def fetch_sequence(self, identifier: str) -> Optional[SeqRecord]:
        """
        Fetch a sequence from Addgene by plasmid ID
        
        Args:
            identifier: Addgene plasmid ID (numeric)
            
        Returns:
            SeqRecord object or None if not found
        """
        try:
            logger.info(f"Fetching sequence from Addgene: {identifier}")
            
            # Clean the identifier (remove "Addgene-" prefix if present)
            if identifier.lower().startswith("addgene-"):
                identifier = identifier[8:]
            
            # Fetch the sequence
            url = f"{self.api_url}/plasmids/{identifier}/full_sequences/"
            response = requests.get(url)
            response.raise_for_status()
            
            data = response.json()
            if not data:
                logger.warning(f"No sequence data found for Addgene ID: {identifier}")
                return None
            
            # Extract sequence data
            sequence = data[0]["sequence"]
            name = data[0]["name"]
            
            # Create SeqRecord
            from Bio.Seq import Seq
            record = SeqRecord(
                Seq(sequence),
                id=f"Addgene-{identifier}",
                name=name,
                description=f"Addgene plasmid #{identifier}: {name}"
            )
            
            # Add annotations
            record.annotations = {
                'molecule_type': 'DNA',
                'topology': 'circular',
                'source': 'Addgene'
            }
            
            logger.info(f"Successfully fetched Addgene-{identifier} ({len(record.seq)} bp)")
            return record
            
        except Exception as e:
            logger.error(f"Failed to fetch Addgene-{identifier}: {str(e)}")
            return None
    
    def search_sequence(self, query: str, max_results: int = 10) -> List[Dict]:
        """
        Search for sequences in Addgene matching the query
        
        Args:
            query: Search term
            max_results: Maximum number of results to return
            
        Returns:
            List of dictionaries with sequence metadata
        """
        try:
            logger.info(f"Searching Addgene for: {query}")
            
            # Perform the search
            url = f"{self.api_url}/plasmids/search/?q={query}&limit={max_results}"
            response = requests.get(url)
            response.raise_for_status()
            
            data = response.json()
            results = []
            
            for item in data.get("results", []):
                results.append({
                    "id": f"Addgene-{item['id']}",
                    "title": item.get("name", ""),
                    "length": item.get("length", 0),
                    "source": "Addgene",
                    "type": "Plasmid",
                    "description": item.get("description", "")
                })
            
            logger.info(f"Found {len(results)} results for query: {query}")
            return results
            
        except Exception as e:
            logger.error(f"Error searching Addgene: {str(e)}")
            return []

class iGEMRepository(SequenceRepository):
    """Repository for fetching sequences from iGEM Registry"""
    
    def __init__(self):
        self.base_url = "http://parts.igem.org/partsdb/get_part.cgi"
    
    def fetch_sequence(self, identifier: str) -> Optional[SeqRecord]:
        """
        Fetch a sequence from iGEM by part ID
        
        Args:
            identifier: iGEM part ID (e.g., BBa_K123456)
            
        Returns:
            SeqRecord object or None if not found
        """
        try:
            logger.info(f"Fetching sequence from iGEM: {identifier}")
            
            # Clean the identifier
            if not identifier.upper().startswith("BBa_"):
                identifier = f"BBa_{identifier}"
            
            # Fetch the sequence
            url = f"{self.base_url}?part={identifier}&format=fasta"
            response = requests.get(url)
            
            # Check if response contains valid FASTA
            if not response.text.startswith(">"):
                logger.warning(f"No valid sequence found for iGEM part: {identifier}")
                return None
            
            # Parse the FASTA
            record = SeqIO.read(io.StringIO(response.text), "fasta")
            
            # Add annotations
            record.annotations = {
                'molecule_type': 'DNA',
                'source': 'iGEM'
            }
            
            logger.info(f"Successfully fetched {identifier} from iGEM ({len(record.seq)} bp)")
            return record
            
        except Exception as e:
            logger.error(f"Failed to fetch {identifier} from iGEM: {str(e)}")
            return None
    
    def search_sequence(self, query: str, max_results: int = 10) -> List[Dict]:
        """
        Search for sequences in iGEM Registry matching the query
        
        Args:
            query: Search term
            max_results: Maximum number of results to return
            
        Returns:
            List of dictionaries with sequence metadata
        """
        # Note: iGEM doesn't have a clean API, so this is a placeholder
        # In production, you might need to scrape the website or use a different approach
        logger.warning("iGEM search not implemented - would require web scraping")
        return []

class LocalRepository(SequenceRepository):
    """Repository for managing local sequence files"""
    
    def __init__(self, data_dir: str = "data"):
        self.data_dir = data_dir
        self.index = self._build_index()
    
    def _build_index(self) -> Dict[str, str]:
        """Build an index of local sequence files"""
        index = {}
        
        if not os.path.exists(self.data_dir):
            logger.warning(f"Data directory not found: {self.data_dir}")
            return index
        
        # Index all sequence files
        for filename in os.listdir(self.data_dir):
            filepath = os.path.join(self.data_dir, filename)
            if os.path.isfile(filepath) and (filename.endswith(".gb") or filename.endswith(".fasta")):
                # Extract ID from filename
                sequence_id = os.path.splitext(filename)[0]
                index[sequence_id.lower()] = filepath
        
        logger.info(f"Indexed {len(index)} local sequence files")
        return index
    
    def fetch_sequence(self, identifier: str) -> Optional[SeqRecord]:
        """
        Fetch a sequence from local files
        
        Args:
            identifier: Sequence identifier
            
        Returns:
            SeqRecord object or None if not found
        """
        try:
            logger.info(f"Fetching sequence from local repository: {identifier}")
            
            # Normalize identifier
            identifier = identifier.lower()
            
            # Check if identifier exists in index
            if identifier not in self.index:
                logger.warning(f"Sequence not found in local repository: {identifier}")
                return None
            
            # Get file path
            filepath = self.index[identifier]
            
            # Determine format from file extension
            file_format = "genbank" if filepath.endswith(".gb") else "fasta"
            
            # Read the sequence
            record = SeqIO.read(filepath, file_format)
            
            # Ensure molecule_type is set
            if not hasattr(record, 'annotations') or not record.annotations:
                record.annotations = {}
            if 'molecule_type' not in record.annotations:
                record.annotations['molecule_type'] = 'DNA'
            
            logger.info(f"Successfully fetched {identifier} from local repository ({len(record.seq)} bp)")
            return record
            
        except Exception as e:
            logger.error(f"Failed to fetch {identifier} from local repository: {str(e)}")
            return None
    
    def search_sequence(self, query: str, max_results: int = 10) -> List[Dict]:
        """
        Search for sequences in local repository matching the query
        
        Args:
            query: Search term
            max_results: Maximum number of results to return
            
        Returns:
            List of dictionaries with sequence metadata
        """
        try:
            logger.info(f"Searching local repository for: {query}")
            
            query = query.lower()
            results = []
            
            for id, filepath in self.index.items():
                if query in id:
                    # Get basic info about the sequence
                    file_format = "genbank" if filepath.endswith(".gb") else "fasta"
                    record = SeqIO.read(filepath, file_format)
                    
                    results.append({
                        "id": id,
                        "title": record.name,
                        "length": len(record.seq),
                        "source": "Local",
                        "type": "DNA",
                        "description": record.description
                    })
                    
                    if len(results) >= max_results:
                        break
            
            logger.info(f"Found {len(results)} results for query: {query}")
            return results
            
        except Exception as e:
            logger.error(f"Error searching local repository: {str(e)}")
            return []
    
    def save_sequence(self, record: SeqRecord, formats: List[str] = ["gb", "fasta"]) -> bool:
        """
        Save a sequence to the local repository
        
        Args:
            record: SeqRecord object to save
            formats: List of formats to save (gb, fasta)
            
        Returns:
            True if successful, False otherwise
        """
        try:
            logger.info(f"Saving sequence to local repository: {record.id}")
            
            # Ensure data directory exists
            os.makedirs(self.data_dir, exist_ok=True)
            
            # Save in each requested format
            for fmt in formats:
                if fmt == "gb" or fmt == "genbank":
                    filepath = os.path.join(self.data_dir, f"{record.id}.gb")
                    SeqIO.write(record, filepath, "genbank")
                elif fmt == "fasta":
                    filepath = os.path.join(self.data_dir, f"{record.id}.fasta")
                    SeqIO.write(record, filepath, "fasta")
            
            # Update index
            self.index[record.id.lower()] = os.path.join(self.data_dir, f"{record.id}.{formats[0]}")
            
            logger.info(f"Successfully saved {record.id} to local repository")
            return True
            
        except Exception as e:
            logger.error(f"Failed to save {record.id} to local repository: {str(e)}")
            return False

class SequenceDatabase:
    """
    Unified sequence database that fetches from multiple repositories
    """
    
    def __init__(self, data_dir: str = "data"):
        self.repositories = {
            "local": LocalRepository(data_dir),
            "ncbi": NCBIRepository(),
            "addgene": AddgeneRepository(),
            "igem": iGEMRepository()
        }
        self.cache = {}
    
    def get_sequence(self, identifier: str, repository: str = None) -> Optional[SeqRecord]:
        """
        Get a sequence by identifier, optionally from a specific repository
        
        Args:
            identifier: Sequence identifier
            repository: Optional repository name to query specifically
            
        Returns:
            SeqRecord object or None if not found
        """
        # Check cache first
        cache_key = f"{repository or 'any'}:{identifier}"
        if cache_key in self.cache:
            logger.info(f"Cache hit for {cache_key}")
            return self.cache[cache_key]
        
        # If repository is specified, query only that one
        if repository and repository in self.repositories:
            record = self.repositories[repository].fetch_sequence(identifier)
            if record:
                self.cache[cache_key] = record
                # Also save to local repository for future use
                self.repositories["local"].save_sequence(record)
                return record
            return None
        
        # Try local repository first
        record = self.repositories["local"].fetch_sequence(identifier)
        if record:
            self.cache[cache_key] = record
            return record
        
        # Try to detect repository from identifier format
        if identifier.startswith("NC_") or identifier.startswith("NM_") or identifier.startswith("XM_"):
            # NCBI accession
            record = self.repositories["ncbi"].fetch_sequence(identifier)
        elif identifier.lower().startswith("addgene-") or identifier.isdigit():
            # Addgene ID
            record = self.repositories["addgene"].fetch_sequence(identifier)
        elif identifier.upper().startswith("BBA_") or (len(identifier) > 2 and identifier[0] == 'B' and identifier[1] == 'B'):
            # iGEM part
            record = self.repositories["igem"].fetch_sequence(identifier)
        else:
            # Try all repositories in order
            for repo_name, repo in self.repositories.items():
                if repo_name == "local":
                    continue  # Already tried
                record = repo.fetch_sequence(identifier)
                if record:
                    break
        
        # Cache and save the result if found
        if record:
            self.cache[cache_key] = record
            self.repositories["local"].save_sequence(record)
        
        return record
    
    def search_sequences(self, query: str, repositories: List[str] = None, max_results: int = 10) -> List[Dict]:
        """
        Search for sequences across repositories
        
        Args:
            query: Search term
            repositories: List of repositories to search (defaults to all)
            max_results: Maximum results per repository
            
        Returns:
            List of dictionaries with sequence metadata
        """
        results = []
        
        # Determine which repositories to search
        repos_to_search = repositories or self.repositories.keys()
        
        # Search each repository
        for repo_name in repos_to_search:
            if repo_name in self.repositories:
                repo_results = self.repositories[repo_name].search_sequence(query, max_results)
                results.extend(repo_results)
        
        return results
    
    def get_available_sequences(self, sequence_type: str = None) -> List[Dict]:
        """
        Get a list of all available sequences in the local repository
        
        Args:
            sequence_type: Optional filter by sequence type
            
        Returns:
            List of dictionaries with sequence metadata
        """
        return self.repositories["local"].search_sequence("", max_results=1000)

    def _fetch_from_ncbi(self, accession):
        """
        Fetch a sequence from NCBI
        
        Args:
            accession: NCBI accession number
            
        Returns:
            SeqRecord object if found, None otherwise
        """
        try:
            # Search for the accession
            logger.info(f"Searching NCBI for accession {accession}")
            handle = Entrez.esearch(db="nucleotide", term=accession)
            record = Entrez.read(handle)
            handle.close()
            
            if not record["IdList"]:
                logger.warning(f"No NCBI records found for {accession}")
                return None
            
            # Get the first result
            id_list = record["IdList"]
            handle = Entrez.efetch(db="nucleotide", id=id_list[0], rettype="gb", retmode="text")
            sequence = SeqIO.read(handle, "genbank")
            handle.close()
            
            logger.info(f"Successfully retrieved {accession} from NCBI")
            
            # Validate sequence
            if self._validate_sequence(sequence):
                # Save locally for future use
                filepath = os.path.join(self.data_dir, f"{accession}.gb")
                SeqIO.write(sequence, filepath, "genbank")
                return sequence
            else:
                logger.warning(f"Sequence {accession} from NCBI failed validation")
                return None
                
        except Exception as e:
            logger.error(f"Error fetching {accession} from NCBI: {str(e)}")
            return None
    
    def _fetch_from_addgene(self, plasmid_id):
        """
        Fetch a sequence from Addgene
        
        Args:
            plasmid_id: Addgene plasmid ID
            
        Returns:
            SeqRecord object if found, None otherwise
        """
        try:
            # Extract numeric ID if needed
            numeric_id = re.search(r'\d+', str(plasmid_id))
            if numeric_id:
                numeric_id = numeric_id.group(0)
            else:
                numeric_id = plasmid_id
            
            # Construct URL for Addgene API
            url = f"https://www.addgene.org/api/v1/plasmids/{numeric_id}/"
            
            logger.info(f"Fetching plasmid {numeric_id} from Addgene")
            
            # Send request
            response = requests.get(url, headers={"Accept": "application/json"})
            
            if response.status_code != 200:
                logger.warning(f"Failed to fetch from Addgene: Status {response.status_code}")
                return None
            
            # Parse JSON response
            data = response.json()
            
            if "sequence" not in data:
                logger.warning("No sequence data in Addgene response")
                return None
            
            # Create SeqRecord from response
            sequence = data["sequence"]
            name = data.get("name", f"Addgene_{numeric_id}")
            description = data.get("description", "")
            
            record = SeqRecord(
                Seq(sequence),
                id=f"Addgene_{numeric_id}",
                name=name,
                description=description
            )
            
            # Add features if available
            if "features" in data:
                self._add_features_from_addgene(record, data["features"])
            
            logger.info(f"Successfully retrieved plasmid {numeric_id} from Addgene")
            
            # Validate sequence
            if self._validate_sequence(record):
                # Save locally for future use
                filepath = os.path.join(self.data_dir, f"Addgene_{numeric_id}.gb")
                SeqIO.write(record, filepath, "genbank")
                return record
            else:
                logger.warning(f"Sequence Addgene_{numeric_id} failed validation")
                return None
                
        except Exception as e:
            logger.error(f"Error fetching plasmid {plasmid_id} from Addgene: {str(e)}")
            return None
    
    def _add_features_from_addgene(self, record, features_data):
        """Add features from Addgene data to SeqRecord"""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        
        for feature in features_data:
            feature_type = feature.get("type", "misc_feature")
            start = feature.get("start", 0)
            end = feature.get("end", 0)
            strand = 1 if feature.get("strand", "+") == "+" else -1
            
            qualifiers = {}
            
            # Add label if available
            if "label" in feature:
                qualifiers["label"] = [feature["label"]]
            
            # Add note if available
            if "note" in feature:
                qualifiers["note"] = [feature["note"]]
            
            # Create feature
            seq_feature = SeqFeature(
                FeatureLocation(start, end),
                type=feature_type,
                qualifiers=qualifiers,
                strand=strand
            )
            
            record.features.append(seq_feature)
    
    def _is_ncbi_accession(self, accession):
        """Check if a string looks like an NCBI accession"""
        # Common NCBI accession patterns
        patterns = [
            r'^[A-Z]{1,2}\d{5,8}$',          # GenBank
            r'^[A-Z]{2}_\d{6,9}$',           # RefSeq
            r'^[A-Z]{3}\d{5}$',              # EMBL/DDBJ
            r'^[A-Z]{4}\d{8,10}$',           # Short reads archive
            r'^[A-Z]{1}\d{5}\.[1-9]$',       # Versioned accession
            r'^NC_\d{6}$',                   # RefSeq chromosome
            r'^NM_\d{6}$'                    # RefSeq mRNA
        ]
        
        for pattern in patterns:
            if re.match(pattern, str(accession)):
                return True
                
        return False
    
    def _is_addgene_id(self, plasmid_id):
        """Check if a string looks like an Addgene ID"""
        # Check for Addgene format (e.g., "Addgene_12345" or just a number)
        if re.match(r'^Addgene_\d+$', str(plasmid_id)):
            return True
        
        # Check if it's a numeric ID between typical Addgene ranges
        if re.match(r'^\d+$', str(plasmid_id)):
            numeric_id = int(plasmid_id)
            return 1 <= numeric_id <= 200000  # Typical Addgene ID range
        
        return False
    
    def _validate_sequence(self, record):
        """
        Validate a sequence record
        
        Args:
            record: SeqRecord object
            
        Returns:
            Boolean indicating if sequence is valid
        """
        # Check if record is valid
        if not isinstance(record, SeqRecord):
            logger.warning("Not a valid SeqRecord object")
            return False
        
        # Check if sequence exists
        if not hasattr(record, 'seq') or not record.seq:
            logger.warning(f"Record {record.id} has no sequence")
            return False
        
        # Check if sequence has valid length
        if len(record.seq) == 0:
            logger.warning(f"Record {record.id} has zero length")
            return False
        
        # Check for excessive ambiguous bases
        ambiguous_count = str(record.seq).count('N') + str(record.seq).count('X')
        if ambiguous_count > len(record.seq) * 0.1:  # More than 10% ambiguous
            logger.warning(f"Record {record.id} has too many ambiguous bases ({ambiguous_count})")
            return False
        
        # Check for valid DNA alphabet
        invalid_chars = set(str(record.seq).upper()) - set('ATGCNRYKMSWBDHV')
        if invalid_chars:
            logger.warning(f"Record {record.id} contains invalid characters: {invalid_chars}")
            return False
        
        # Add more checks as needed
        
        return True
    
    def create_record_from_sequence(self, sequence, sequence_id, name=None, description=None):
        """
        Create a SeqRecord from a raw sequence string
        
        Args:
            sequence: Sequence string
            sequence_id: Identifier for the sequence
            name: Optional name
            description: Optional description
            
        Returns:
            SeqRecord object if successful, None otherwise
        """
        try:
            # Validate sequence string
            sequence = sequence.strip().upper()
            invalid_chars = set(sequence) - set('ATGCNRYKMSWBDHV')
            
            if invalid_chars:
                logger.warning(f"Sequence contains invalid characters: {invalid_chars}")
                return None
            
            # Create SeqRecord
            record = SeqRecord(
                Seq(sequence),
                id=sequence_id,
                name=name or sequence_id,
                description=description or ""
            )
            
            # Add annotations
            record.annotations["molecule_type"] = "DNA"
            
            # Validate
            if self._validate_sequence(record):
                return record
            else:
                return None
                
        except Exception as e:
            logger.error(f"Error creating record: {str(e)}")
            return None
    
    def import_sequence_from_string(self, sequence_string, format_type="fasta"):
        """
        Import a sequence from a formatted string
        
        Args:
            sequence_string: String containing sequence in specified format
            format_type: Format of the sequence (fasta, genbank)
            
        Returns:
            Imported SeqRecord object if successful, None otherwise
        """
        try:
            # Parse from string
            handle = StringIO(sequence_string)
            record = SeqIO.read(handle, format_type)
            
            # Validate
            if self._validate_sequence(record):
                # Add to database
                self.sequences[record.id] = record
                return record
            else:
                logger.warning(f"Imported sequence failed validation")
                return None
                
        except Exception as e:
            logger.error(f"Error importing sequence: {str(e)}")
            return None
    
    def search_sequences(self, query, search_description=True):
        """
        Search for sequences matching a query
        
        Args:
            query: Search query
            search_description: Whether to search in descriptions
            
        Returns:
            List of matching sequence IDs
        """
        query = query.lower()
        results = []
        
        for seq_id, record in self.sequences.items():
            # Search in ID and name
            if query in seq_id.lower() or (hasattr(record, 'name') and query in record.name.lower()):
                results.append(seq_id)
                continue
            
            # Search in description
            if search_description and hasattr(record, 'description') and query in record.description.lower():
                results.append(seq_id)
                continue
            
            # Search in features
            if hasattr(record, 'features'):
                for feature in record.features:
                    if not hasattr(feature, 'qualifiers'):
                        continue
                        
                    for qualifier_values in feature.qualifiers.values():
                        for value in qualifier_values:
                            if query in value.lower():
                                results.append(seq_id)
                                break
        
        return results

# Initialize the database
sequence_db = SequenceDatabase()

def get_sequence(identifier: str, repository: str = None) -> Optional[SeqRecord]:
    """
    Get a sequence by identifier (convenience function)
    
    Args:
        identifier: Sequence identifier
        repository: Optional repository name
        
    Returns:
        SeqRecord object or None if not found
    """
    return sequence_db.get_sequence(identifier, repository)

def search_sequences(query: str, repositories: List[str] = None, max_results: int = 10) -> List[Dict]:
    """
    Search for sequences (convenience function)
    
    Args:
        query: Search term
        repositories: List of repositories to search
        max_results: Maximum results per repository
        
    Returns:
        List of dictionaries with sequence metadata
    """
    return sequence_db.search_sequences(query, repositories, max_results)

def get_available_sequences(sequence_type: str = None) -> List[Dict]:
    """
    Get available sequences (convenience function)
    
    Args:
        sequence_type: Optional filter by sequence type
        
    Returns:
        List of dictionaries with sequence metadata
    """
    return sequence_db.get_available_sequences(sequence_type) 
# Privacy and Data Handling

This Shiny application processes user-provided count and metadata files within the active session.  
Uploaded files are handled in temporary storage and are not retained by the application after the session ends.

If optional pathway analysis is enabled, the derived list of differentially expressed gene symbols is submitted to the external Enrichr service for enrichment analysis. The uploaded count matrix and metadata files themselves are not sent to Enrichr.

### Deployment context

The app is hosted at **https://r-odaf.nl** inside a Docker container managed by the R-ODAF team at Maastricht University.  
Standard server access logs, such as IP address, timestamp, and browser information, may be retained temporarily for technical monitoring.

### User responsibility

Users should avoid uploading personally identifiable or confidential data.  
The application is intended for research RNA-Seq count matrices and associated sample metadata.